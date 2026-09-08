use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};

use crate::store::{GeneInterner, Packer, TxMeta};

pub const SHARD_COUNT: u32 = 16;

#[derive(Debug, thiserror::Error)]
pub enum FastaError {
    #[error("missing shard {0}")]
    MissingShard(PathBuf),
    #[error("io error reading {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("{0}")]
    Other(String),
}

pub fn shard_path(data_dir: &Path, i: u32) -> PathBuf {
    data_dir.join(format!("human.{i}.rna.fna"))
}

pub fn required_shards(data_dir: &Path) -> Result<Vec<PathBuf>, FastaError> {
    let mut out = Vec::with_capacity(SHARD_COUNT as usize);
    for i in 1..=SHARD_COUNT {
        let p = shard_path(data_dir, i);
        if !p.is_file() {
            return Err(FastaError::MissingShard(p));
        }
        out.push(p);
    }
    Ok(out)
}

/// Fast fingerprint: name + size + mtime for each required shard.
pub fn fasta_fingerprint(paths: &[PathBuf]) -> Result<String, FastaError> {
    let mut hasher = Sha256::new();
    for p in paths {
        let meta = std::fs::metadata(p).map_err(|source| FastaError::Io {
            path: p.clone(),
            source,
        })?;
        hasher.update(
            p.file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .as_bytes(),
        );
        hasher.update(b"\0");
        hasher.update(meta.len().to_le_bytes());
        let mtime = meta
            .modified()
            .ok()
            .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
            .map(|d| d.as_secs())
            .unwrap_or(0);
        hasher.update(mtime.to_le_bytes());
    }
    Ok(hex::encode(hasher.finalize()))
}

/// Parse `>ACCESSION ... (GENE), ..., MOL_TYPE`.
pub fn parse_header(line: &str) -> (String, String, String) {
    let line = line.strip_prefix('>').unwrap_or(line).trim();
    let accession = line.split_whitespace().next().unwrap_or("").to_string();
    let mol_type = line
        .rsplit_once(',')
        .map(|(_, t)| t.trim().to_string())
        .unwrap_or_default();
    let mut gene = String::new();
    let mut start = None;
    for (i, c) in line.char_indices() {
        if c == '(' {
            start = Some(i + 1);
        } else if c == ')' {
            if let Some(s) = start {
                gene = line[s..i].to_string();
            }
        }
    }
    (accession, gene, mol_type)
}

pub fn strip_version(accession: &str) -> &str {
    match accession.rsplit_once('.') {
        Some((stem, ver)) if ver.bytes().all(|b| b.is_ascii_digit()) => stem,
        _ => accession,
    }
}

pub fn load_fastas(
    paths: &[PathBuf],
    packer: &mut Packer,
    metas: &mut Vec<TxMeta>,
    genes: &mut GeneInterner,
) -> Result<(), FastaError> {
    for path in paths {
        load_one(path, packer, metas, genes)?;
    }
    packer.pad();
    Ok(())
}

fn load_one(
    path: &Path,
    packer: &mut Packer,
    metas: &mut Vec<TxMeta>,
    genes: &mut GeneInterner,
) -> Result<(), FastaError> {
    let file = File::open(path).map_err(|source| FastaError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let reader = BufReader::with_capacity(1 << 20, file);
    let mut header: Option<String> = None;
    let mut seq = String::new();

    for line in reader.lines() {
        let line = line.map_err(|source| FastaError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(h) = header.take() {
                flush_record(&h, &seq, packer, metas, genes);
                seq.clear();
            }
            header = Some(rest.to_string());
        } else {
            seq.push_str(line.trim());
        }
    }
    if let Some(h) = header {
        flush_record(&h, &seq, packer, metas, genes);
    }
    Ok(())
}

fn flush_record(
    header: &str,
    seq: &str,
    packer: &mut Packer,
    metas: &mut Vec<TxMeta>,
    genes: &mut GeneInterner,
) {
    let (accession, gene, mol_type) = parse_header(header);
    if accession.is_empty() || seq.is_empty() {
        return;
    }
    let gene_key = gene.to_ascii_uppercase();
    let gene_id = genes.intern(&gene_key);
    let accession_key = strip_version(&accession).to_ascii_uppercase();
    let (offset, len, has_n) = packer.push_seq(seq.as_bytes());
    metas.push(TxMeta {
        accession,
        accession_key,
        gene,
        gene_key,
        mol_type,
        offset,
        len,
        gene_id,
        has_n,
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_pcsk9_style() {
        let h =
            ">NM_174936.4 Homo sapiens proprotein convertase subtilisin/kexin type 9 (PCSK9), mRNA";
        let (acc, gene, mol) = parse_header(h);
        assert_eq!(acc, "NM_174936.4");
        assert_eq!(gene, "PCSK9");
        assert_eq!(mol, "mRNA");
    }

    #[test]
    fn header_variant() {
        let h = ">NM_001438777.1 Homo sapiens caldesmon 1 (CALD1), transcript variant 18, mRNA";
        let (acc, gene, mol) = parse_header(h);
        assert_eq!(acc, "NM_001438777.1");
        assert_eq!(gene, "CALD1");
        assert_eq!(mol, "mRNA");
    }

    #[test]
    fn strip_acc_version() {
        assert_eq!(strip_version("NM_174936.4"), "NM_174936");
        assert_eq!(strip_version("NM_174936"), "NM_174936");
    }
}
