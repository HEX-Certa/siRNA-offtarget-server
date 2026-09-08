use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use tracing::{info, warn};

use crate::fasta::{self, FastaError};
use crate::seq::encode_base;

const SNAPSHOT_MAGIC: &[u8; 8] = b"SOTCv1\n\0";

#[derive(Debug, thiserror::Error)]
pub enum StoreError {
    #[error(transparent)]
    Fasta(#[from] FastaError),
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("invalid snapshot: {0}")]
    Snapshot(String),
}

#[derive(Debug, Clone)]
pub struct TxMeta {
    pub accession: String,
    pub accession_key: String,
    pub gene: String,
    pub gene_key: String,
    pub mol_type: String,
    pub offset: u64,
    pub len: u32,
    pub gene_id: u32,
    pub has_n: bool,
}

#[derive(Debug, Default)]
pub struct GeneInterner {
    map: HashMap<String, u32>,
    names: Vec<String>,
}

impl GeneInterner {
    pub fn intern(&mut self, key: &str) -> u32 {
        if let Some(&id) = self.map.get(key) {
            return id;
        }
        let id = self.names.len() as u32;
        self.map.insert(key.to_string(), id);
        self.names.push(key.to_string());
        id
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }
}

#[derive(Debug, Default)]
pub struct Packer {
    pub packed: Vec<u64>,
    pub nmask: Vec<u64>,
    pub bases: u64,
}

impl Packer {
    pub fn push_seq(&mut self, seq: &[u8]) -> (u64, u32, bool) {
        let offset = self.bases;
        let mut has_n = false;
        for &b in seq {
            if b.is_ascii_whitespace() {
                continue;
            }
            let (bits, is_n) = encode_base(b);
            has_n |= is_n;
            self.push_base(bits, is_n);
        }
        let len = (self.bases - offset) as u32;
        (offset, len, has_n)
    }

    fn push_base(&mut self, bits: u8, is_n: bool) {
        let i = self.bases;
        let w = (i / 32) as usize;
        let shift = ((i % 32) * 2) as u32;
        if w >= self.packed.len() {
            self.packed.push(0);
        }
        self.packed[w] |= (bits as u64) << shift;

        let nw = (i / 64) as usize;
        let nshift = (i % 64) as u32;
        if nw >= self.nmask.len() {
            self.nmask.push(0);
        }
        if is_n {
            self.nmask[nw] |= 1u64 << nshift;
        }
        self.bases += 1;
    }

    pub fn pad(&mut self) {
        self.packed.push(0);
        self.nmask.push(0);
    }
}

#[derive(Debug)]
pub struct Store {
    pub packed: Vec<u64>,
    pub nmask: Vec<u64>,
    pub metas: Vec<TxMeta>,
    pub gene_names: Vec<String>,
    pub bases: u64,
    pub fingerprint: String,
    pub shards: Vec<String>,
    /// Snapshot file mtime (UNIX seconds), if known.
    pub indexed_at: Option<u64>,
}

impl Store {
    pub fn transcripts(&self) -> u64 {
        self.metas.len() as u64
    }

    pub fn includes_xm_xr(&self) -> bool {
        self.metas.iter().any(|tx| {
            let a = tx.accession.as_str();
            a.starts_with("XM_") || a.starts_with("XR_")
        })
    }

    pub fn from_fastas(paths: &[PathBuf], fingerprint: String) -> Result<Self, StoreError> {
        let mut packer = Packer::default();
        let mut metas = Vec::new();
        let mut genes = GeneInterner::default();
        info!(files = paths.len(), "parsing FASTA shards");
        fasta::load_fastas(paths, &mut packer, &mut metas, &mut genes)?;
        let shards = paths
            .iter()
            .map(|p| {
                p.file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("")
                    .to_string()
            })
            .collect();
        info!(
            transcripts = metas.len(),
            bases = packer.bases,
            "FASTA loaded"
        );
        Ok(Self {
            packed: packer.packed,
            nmask: packer.nmask,
            metas,
            gene_names: genes.names().to_vec(),
            bases: packer.bases,
            fingerprint,
            shards,
            indexed_at: None,
        })
    }

    pub fn load_or_build(data_dir: &Path, index_dir: &Path) -> Result<Self, StoreError> {
        let paths = fasta::required_shards(data_dir)?;
        let fp = fasta::fasta_fingerprint(&paths)?;
        fs::create_dir_all(index_dir)?;
        let snap = snapshot_path(index_dir);
        if snap.is_file() {
            match Self::read_snapshot(&snap) {
                Ok(mut store) if store.fingerprint == fp => {
                    store.indexed_at = file_mtime_secs(&snap);
                    info!(path = %snap.display(), "loaded transcriptome snapshot");
                    return Ok(store);
                }
                Ok(store) => {
                    warn!(
                        old = %store.fingerprint,
                        new = %fp,
                        "snapshot fingerprint mismatch, rebuilding"
                    );
                }
                Err(e) => warn!(error = %e, "snapshot unreadable, rebuilding"),
            }
        }
        let mut store = Self::from_fastas(&paths, fp)?;
        store.write_snapshot(&snap)?;
        store.indexed_at = file_mtime_secs(&snap);
        info!(path = %snap.display(), "wrote transcriptome snapshot");
        Ok(store)
    }

    pub fn write_snapshot(&self, path: &Path) -> Result<(), StoreError> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        let mut f = File::create(path)?;
        f.write_all(SNAPSHOT_MAGIC)?;
        write_u64(&mut f, self.metas.len() as u64)?;
        write_u64(&mut f, self.bases)?;
        write_u64(&mut f, self.packed.len() as u64)?;
        write_u64(&mut f, self.nmask.len() as u64)?;
        write_str(&mut f, &self.fingerprint)?;
        write_u32(&mut f, self.shards.len() as u32)?;
        for s in &self.shards {
            write_str(&mut f, s)?;
        }
        write_u32(&mut f, self.gene_names.len() as u32)?;
        for s in &self.gene_names {
            write_str(&mut f, s)?;
        }
        for w in &self.packed {
            f.write_all(&w.to_le_bytes())?;
        }
        for w in &self.nmask {
            f.write_all(&w.to_le_bytes())?;
        }
        for tx in &self.metas {
            write_str(&mut f, &tx.accession)?;
            write_str(&mut f, &tx.accession_key)?;
            write_str(&mut f, &tx.gene)?;
            write_str(&mut f, &tx.gene_key)?;
            write_str(&mut f, &tx.mol_type)?;
            write_u64(&mut f, tx.offset)?;
            write_u32(&mut f, tx.len)?;
            write_u32(&mut f, tx.gene_id)?;
            f.write_all(&[u8::from(tx.has_n)])?;
        }
        f.flush()?;
        Ok(())
    }

    pub fn read_snapshot(path: &Path) -> Result<Self, StoreError> {
        let mut f = File::open(path)?;
        let mut magic = [0u8; 8];
        f.read_exact(&mut magic)?;
        if &magic != SNAPSHOT_MAGIC {
            return Err(StoreError::Snapshot("bad magic".into()));
        }
        let n_tx = read_u64(&mut f)? as usize;
        let bases = read_u64(&mut f)?;
        let n_packed = read_u64(&mut f)? as usize;
        let n_nmask = read_u64(&mut f)? as usize;
        let fingerprint = read_str(&mut f)?;
        let n_shards = read_u32(&mut f)? as usize;
        let mut shards = Vec::with_capacity(n_shards);
        for _ in 0..n_shards {
            shards.push(read_str(&mut f)?);
        }
        let n_genes = read_u32(&mut f)? as usize;
        let mut gene_names = Vec::with_capacity(n_genes);
        for _ in 0..n_genes {
            gene_names.push(read_str(&mut f)?);
        }
        let mut packed = vec![0u64; n_packed];
        for w in &mut packed {
            *w = read_u64(&mut f)?;
        }
        let mut nmask = vec![0u64; n_nmask];
        for w in &mut nmask {
            *w = read_u64(&mut f)?;
        }
        let mut metas = Vec::with_capacity(n_tx);
        for _ in 0..n_tx {
            metas.push(TxMeta {
                accession: read_str(&mut f)?,
                accession_key: read_str(&mut f)?,
                gene: read_str(&mut f)?,
                gene_key: read_str(&mut f)?,
                mol_type: read_str(&mut f)?,
                offset: read_u64(&mut f)?,
                len: read_u32(&mut f)?,
                gene_id: read_u32(&mut f)?,
                has_n: {
                    let mut b = [0u8; 1];
                    f.read_exact(&mut b)?;
                    b[0] != 0
                },
            });
        }
        Ok(Self {
            packed,
            nmask,
            metas,
            gene_names,
            bases,
            fingerprint,
            shards,
            indexed_at: None,
        })
    }
}

pub fn snapshot_path(index_dir: &Path) -> PathBuf {
    index_dir.join("transcriptome.bin")
}

fn file_mtime_secs(path: &Path) -> Option<u64> {
    std::fs::metadata(path)
        .ok()
        .and_then(|m| m.modified().ok())
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|d| d.as_secs())
}

fn write_u32(w: &mut impl Write, v: u32) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn write_u64(w: &mut impl Write, v: u64) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn write_str(w: &mut impl Write, s: &str) -> std::io::Result<()> {
    write_u32(w, s.len() as u32)?;
    w.write_all(s.as_bytes())
}

fn read_u32(r: &mut impl Read) -> std::io::Result<u32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64(r: &mut impl Read) -> std::io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_str(r: &mut impl Read) -> Result<String, StoreError> {
    let n = read_u32(r)? as usize;
    if n > 1_000_000 {
        return Err(StoreError::Snapshot("string too long".into()));
    }
    let mut buf = vec![0u8; n];
    r.read_exact(&mut buf)?;
    String::from_utf8(buf).map_err(|e| StoreError::Snapshot(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::decode_window;

    #[test]
    #[ignore = "loads the full RefSeq transcriptome"]
    fn smoke_full_transcriptome_random_21mer() {
        let data = std::path::Path::new("data");
        let tmp = tempfile::tempdir().unwrap();
        let store = Store::load_or_build(data, tmp.path()).expect("load data/");
        assert!(store.transcripts() > 200_000);
        let items = [crate::search::SearchItem {
            id: "rand".into(),
            guide: "GATCGTAGCTAGGCTTAGCTA".into(),
            sense: None,
        }];
        let r = crate::search::search_batch(
            &store,
            &items,
            &crate::search::TargetFilter::new(None, &[]),
            crate::search::ScanMode::T6b { scan_sense: false },
        );
        assert_eq!(
            r[0].offtarget.as_ref().unwrap().perfect,
            0,
            "random 21-mer should have no perfect transcriptome hit"
        );
    }

    #[test]
    fn pack_and_snapshot_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("tiny.fna");
        std::fs::write(
            &fa,
            ">NM_000001.1 Homo sapiens dummy (DUMMY), mRNA\nACGTACGTAAATAAACNCCAGG\n",
        )
        .unwrap();
        let store = Store::from_fastas(&[fa], "fp1".into()).unwrap();
        assert_eq!(store.transcripts(), 1);
        assert_eq!(store.metas[0].gene, "DUMMY");
        assert!(store.metas[0].has_n);
        assert_eq!(decode_window(&store.packed, 0, 4), "ACGT");

        let snap = dir.path().join("snap.bin");
        store.write_snapshot(&snap).unwrap();
        let loaded = Store::read_snapshot(&snap).unwrap();
        assert_eq!(loaded.fingerprint, "fp1");
        assert_eq!(loaded.metas[0].accession, "NM_000001.1");
        assert_eq!(loaded.bases, store.bases);
        assert_eq!(loaded.packed, store.packed);
    }
}
