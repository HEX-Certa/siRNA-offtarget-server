use std::collections::{HashMap, HashSet};

use rayon::prelude::*;

use crate::fasta::strip_version;
use crate::score::specificity_label;
use crate::seq::{
    decode_window, extract_packed, hamming_packed, oligo_2_20, pack_oligo, revcomp, seed_of_guide,
    window_has_n, MAX_MISMATCH, SEED_LEN,
};
use crate::store::{Store, TxMeta};
use crate::types::{
    Hit, OfftargetCounts, OnTargetCounts, Query19mer, QueryResult, StrandMinMismatch, StrandMmBins,
    StrandResult, INTERNAL_HIT_CAP, MIN_MISMATCH_NONE,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScanMode {
    T6b { scan_sense: bool },
    Sidirect,
}

#[derive(Debug, Clone)]
pub struct TargetFilter {
    gene: Option<String>,
    accessions: HashSet<String>,
}

impl TargetFilter {
    pub fn new(gene: Option<&str>, accessions: &[String]) -> Self {
        let gene = gene
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|s| s.to_ascii_uppercase());
        let mut set = HashSet::new();
        for a in accessions {
            let n = a.trim().to_ascii_uppercase();
            if n.is_empty() {
                continue;
            }
            set.insert(strip_version(&n).to_string());
            set.insert(n);
        }
        Self {
            gene,
            accessions: set,
        }
    }

    pub fn is_on_target(&self, tx: &TxMeta) -> bool {
        if let Some(g) = &self.gene {
            if tx.gene_key == *g {
                return true;
            }
        }
        if self.accessions.contains(&tx.accession.to_ascii_uppercase())
            || self.accessions.contains(&tx.accession_key)
        {
            return true;
        }
        false
    }
}

#[derive(Clone)]
struct OligoPattern {
    site: u64,
    seed: u64,
    len: usize,
    scan_seed: bool,
}

impl OligoPattern {
    fn from_full_guide(guide: &str) -> Self {
        let site_seq = revcomp(guide);
        let seed_seq = revcomp(seed_of_guide(guide).expect("validated guide"));
        Self {
            site: pack_oligo(&site_seq),
            seed: pack_oligo(&seed_seq),
            len: guide.len(),
            scan_seed: true,
        }
    }

    fn from_19mer(oligo19: &str) -> Self {
        let site_seq = revcomp(oligo19);
        Self {
            site: pack_oligo(&site_seq),
            seed: 0,
            len: oligo19.len(),
            scan_seed: false,
        }
    }
}

#[derive(Default)]
struct LocalStrand {
    mm: [u64; 4],
    on_tx: HashSet<u32>,
    on_genes: HashSet<u32>,
    seed_tx: HashSet<u32>,
    seed_genes: HashSet<u32>,
    hits: Vec<Hit>,
}

impl LocalStrand {
    fn add_near(&mut self, hit: Hit, on_target: bool, tx_idx: u32, gene_id: u32) {
        if on_target {
            self.on_tx.insert(tx_idx);
            self.on_genes.insert(gene_id);
        } else if (hit.mismatches as usize) < self.mm.len() {
            self.mm[hit.mismatches as usize] += 1;
        }
        push_hit(&mut self.hits, hit);
    }

    fn add_seed(&mut self, on_target: bool, tx_idx: u32, gene_id: u32, gene_nonempty: bool) {
        if on_target {
            return;
        }
        self.seed_tx.insert(tx_idx);
        if gene_nonempty {
            self.seed_genes.insert(gene_id);
        }
    }

    fn merge(&mut self, other: LocalStrand) {
        for i in 0..4 {
            self.mm[i] += other.mm[i];
        }
        self.on_tx.extend(other.on_tx);
        self.on_genes.extend(other.on_genes);
        self.seed_tx.extend(other.seed_tx);
        self.seed_genes.extend(other.seed_genes);
        for h in other.hits {
            push_hit(&mut self.hits, h);
        }
    }

    fn finish(mut self) -> (OnTargetCounts, OfftargetCounts, Vec<Hit>) {
        self.hits.sort_by(|a, b| {
            a.mismatches
                .cmp(&b.mismatches)
                .then_with(|| a.accession.cmp(&b.accession))
                .then_with(|| a.position.cmp(&b.position))
        });
        if self.hits.len() > INTERNAL_HIT_CAP {
            self.hits.truncate(INTERNAL_HIT_CAP);
        }
        let on_target = OnTargetCounts {
            transcripts: self.on_tx.len() as u64,
            genes: self.on_genes.len() as u64,
        };
        let offtarget = OfftargetCounts {
            perfect: self.mm[0],
            mismatch_1: self.mm[1],
            mismatch_2: self.mm[2],
            mismatch_3: self.mm[3],
            seed_transcripts: self.seed_tx.len() as u64,
            seed_genes: self.seed_genes.len() as u64,
        };
        (on_target, offtarget, self.hits)
    }
}

/// siDirect: one count per transcript (min Hamming), plus gene-level mins.
struct LocalSidirectStrand {
    off_mm: [u64; 4],
    all_mm: [u64; 4],
    on_tx: HashSet<u32>,
    on_genes: HashSet<u32>,
    gene_min_off: HashMap<u32, u8>,
    gene_min_all: HashMap<u32, u8>,
    off_min: u8,
    hits: Vec<Hit>,
}

impl LocalSidirectStrand {
    fn new() -> Self {
        Self {
            off_mm: [0; 4],
            all_mm: [0; 4],
            on_tx: HashSet::new(),
            on_genes: HashSet::new(),
            gene_min_off: HashMap::new(),
            gene_min_all: HashMap::new(),
            off_min: MIN_MISMATCH_NONE,
            hits: Vec::new(),
        }
    }

    fn add_tx(
        &mut self,
        min_mm: u8,
        on_target: bool,
        tx_idx: u32,
        gene_id: u32,
        gene_nonempty: bool,
    ) {
        if min_mm > MAX_MISMATCH {
            return;
        }
        let i = min_mm as usize;
        self.all_mm[i] += 1;
        if gene_nonempty {
            self.gene_min_all
                .entry(gene_id)
                .and_modify(|m| *m = (*m).min(min_mm))
                .or_insert(min_mm);
        }
        if on_target {
            self.on_tx.insert(tx_idx);
            self.on_genes.insert(gene_id);
        } else {
            self.off_mm[i] += 1;
            self.off_min = self.off_min.min(min_mm);
            if gene_nonempty {
                self.gene_min_off
                    .entry(gene_id)
                    .and_modify(|m| *m = (*m).min(min_mm))
                    .or_insert(min_mm);
            }
        }
    }

    fn merge(&mut self, other: LocalSidirectStrand) {
        for i in 0..4 {
            self.off_mm[i] += other.off_mm[i];
            self.all_mm[i] += other.all_mm[i];
        }
        self.on_tx.extend(other.on_tx);
        self.on_genes.extend(other.on_genes);
        self.off_min = self.off_min.min(other.off_min);
        for (g, mm) in other.gene_min_off {
            self.gene_min_off
                .entry(g)
                .and_modify(|m| *m = (*m).min(mm))
                .or_insert(mm);
        }
        for (g, mm) in other.gene_min_all {
            self.gene_min_all
                .entry(g)
                .and_modify(|m| *m = (*m).min(mm))
                .or_insert(mm);
        }
        for h in other.hits {
            push_hit(&mut self.hits, h);
        }
    }

    fn finish(mut self) -> SidirectStrandOut {
        self.hits.sort_by(|a, b| {
            a.mismatches
                .cmp(&b.mismatches)
                .then_with(|| a.accession.cmp(&b.accession))
                .then_with(|| a.position.cmp(&b.position))
        });
        if self.hits.len() > INTERNAL_HIT_CAP {
            self.hits.truncate(INTERNAL_HIT_CAP);
        }
        let mut off_genes = [0u64; 4];
        let mut all_genes = [0u64; 4];
        for mm in self.gene_min_off.values() {
            if (*mm as usize) < 4 {
                off_genes[*mm as usize] += 1;
            }
        }
        for mm in self.gene_min_all.values() {
            if (*mm as usize) < 4 {
                all_genes[*mm as usize] += 1;
            }
        }
        SidirectStrandOut {
            on_target: OnTargetCounts {
                transcripts: self.on_tx.len() as u64,
                genes: self.on_genes.len() as u64,
            },
            offtarget_only: self.off_mm,
            including_on_target: self.all_mm,
            offtarget_only_genes: off_genes,
            including_on_target_genes: all_genes,
            min_mismatch: self.off_min,
            hits: self.hits,
        }
    }
}

struct SidirectStrandOut {
    on_target: OnTargetCounts,
    offtarget_only: [u64; 4],
    including_on_target: [u64; 4],
    offtarget_only_genes: [u64; 4],
    including_on_target_genes: [u64; 4],
    min_mismatch: u8,
    hits: Vec<Hit>,
}

fn push_hit(hits: &mut Vec<Hit>, hit: Hit) {
    hits.push(hit);
    if hits.len() > INTERNAL_HIT_CAP * 2 {
        hits.sort_by(|a, b| {
            a.mismatches
                .cmp(&b.mismatches)
                .then_with(|| a.accession.cmp(&b.accession))
                .then_with(|| a.position.cmp(&b.position))
        });
        hits.truncate(INTERNAL_HIT_CAP);
    }
}

struct PreparedQuery {
    id: String,
    guide: String,
    guide_pat: OligoPattern,
    sense: Option<(String, OligoPattern)>,
    /// For sidirect: the 19-mers actually searched.
    guide_19: Option<String>,
    passenger_19: Option<String>,
}

#[derive(Clone)]
pub struct SearchItem {
    pub id: String,
    pub guide: String,
    pub sense: Option<String>,
}

pub fn search_batch(
    store: &Store,
    items: &[SearchItem],
    filter: &TargetFilter,
    mode: ScanMode,
) -> Vec<QueryResult> {
    match mode {
        ScanMode::T6b { scan_sense } => search_batch_t6b(store, items, filter, scan_sense),
        ScanMode::Sidirect => search_batch_sidirect(store, items, filter),
    }
}

fn search_batch_t6b(
    store: &Store,
    items: &[SearchItem],
    filter: &TargetFilter,
    scan_sense: bool,
) -> Vec<QueryResult> {
    let prepared: Vec<PreparedQuery> = items
        .iter()
        .map(|it| {
            let sense = if scan_sense {
                let seq = it.sense.clone().unwrap_or_else(|| revcomp(&it.guide));
                let pat = OligoPattern::from_full_guide(&seq);
                Some((seq, pat))
            } else {
                None
            };
            PreparedQuery {
                id: it.id.clone(),
                guide: it.guide.clone(),
                guide_pat: OligoPattern::from_full_guide(&it.guide),
                sense,
                guide_19: None,
                passenger_19: None,
            }
        })
        .collect();

    let n = prepared.len();
    let folded = store
        .metas
        .par_iter()
        .enumerate()
        .fold(
            || {
                (0..n)
                    .map(|_| (LocalStrand::default(), LocalStrand::default()))
                    .collect::<Vec<_>>()
            },
            |mut acc, (tx_idx, tx)| {
                scan_transcript_t6b(store, tx, tx_idx as u32, &prepared, filter, &mut acc);
                acc
            },
        )
        .reduce(
            || {
                (0..n)
                    .map(|_| (LocalStrand::default(), LocalStrand::default()))
                    .collect()
            },
            |mut a, b| {
                for (i, (g, s)) in b.into_iter().enumerate() {
                    a[i].0.merge(g);
                    a[i].1.merge(s);
                }
                a
            },
        );

    prepared
        .into_iter()
        .zip(folded)
        .map(|(pq, (g, s))| {
            let (on_target, offtarget, hits) = g.finish();
            let label = specificity_label(&offtarget).to_string();
            let sense = pq.sense.map(|(seq, _)| {
                let (on_t, off, hits) = s.finish();
                let sl = specificity_label(&off).to_string();
                StrandResult {
                    sequence: seq,
                    specificity_label: sl,
                    on_target: on_t,
                    offtarget: off,
                    hits,
                }
            });
            QueryResult {
                id: pq.id,
                length: pq.guide.len(),
                guide: pq.guide,
                cached: false,
                profile: "t6b".into(),
                specificity_label: Some(label),
                on_target,
                offtarget: Some(offtarget),
                hits,
                sense,
                query_19mer: None,
                min_mismatch: None,
                passes_hide_less_specific: None,
                offtarget_only: None,
                including_on_target: None,
                offtarget_only_genes: None,
                including_on_target_genes: None,
            }
        })
        .collect()
}

fn search_batch_sidirect(
    store: &Store,
    items: &[SearchItem],
    filter: &TargetFilter,
) -> Vec<QueryResult> {
    let prepared: Vec<PreparedQuery> = items
        .iter()
        .map(|it| {
            let g19 = oligo_2_20(&it.guide)
                .expect("validated guide")
                .to_string();
            let passenger = it.sense.clone().unwrap_or_else(|| revcomp(&it.guide));
            let p19 = oligo_2_20(&passenger)
                .expect("validated passenger")
                .to_string();
            PreparedQuery {
                id: it.id.clone(),
                guide: it.guide.clone(),
                guide_pat: OligoPattern::from_19mer(&g19),
                sense: Some((passenger, OligoPattern::from_19mer(&p19))),
                guide_19: Some(g19),
                passenger_19: Some(p19),
            }
        })
        .collect();

    let n = prepared.len();
    let folded = store
        .metas
        .par_iter()
        .enumerate()
        .fold(
            || {
                (0..n)
                    .map(|_| (LocalSidirectStrand::new(), LocalSidirectStrand::new()))
                    .collect::<Vec<_>>()
            },
            |mut acc, (tx_idx, tx)| {
                scan_transcript_sidirect(store, tx, tx_idx as u32, &prepared, filter, &mut acc);
                acc
            },
        )
        .reduce(
            || {
                (0..n)
                    .map(|_| (LocalSidirectStrand::new(), LocalSidirectStrand::new()))
                    .collect()
            },
            |mut a, b| {
                for (i, (g, s)) in b.into_iter().enumerate() {
                    a[i].0.merge(g);
                    a[i].1.merge(s);
                }
                a
            },
        );

    prepared
        .into_iter()
        .zip(folded)
        .map(|(pq, (g, s))| {
            let guide_out = g.finish();
            let pass_out = s.finish();
            let mut hits = guide_out.hits;
            hits.extend(pass_out.hits);
            hits.sort_by(|a, b| {
                a.mismatches
                    .cmp(&b.mismatches)
                    .then_with(|| a.strand.cmp(&b.strand))
                    .then_with(|| a.accession.cmp(&b.accession))
                    .then_with(|| a.position.cmp(&b.position))
            });
            if hits.len() > INTERNAL_HIT_CAP {
                hits.truncate(INTERNAL_HIT_CAP);
            }
            let min_mismatch = StrandMinMismatch {
                guide: guide_out.min_mismatch,
                passenger: pass_out.min_mismatch,
            };
            let passes = min_mismatch.guide >= 2 && min_mismatch.passenger >= 2;
            QueryResult {
                id: pq.id,
                length: pq.guide.len(),
                guide: pq.guide,
                cached: false,
                profile: "sidirect".into(),
                specificity_label: None,
                on_target: guide_out.on_target,
                offtarget: None,
                hits,
                sense: None,
                query_19mer: Some(Query19mer {
                    guide: pq.guide_19.unwrap_or_default(),
                    passenger: pq.passenger_19.unwrap_or_default(),
                }),
                min_mismatch: Some(min_mismatch),
                passes_hide_less_specific: Some(passes),
                offtarget_only: Some(StrandMmBins {
                    guide: guide_out.offtarget_only,
                    passenger: pass_out.offtarget_only,
                }),
                including_on_target: Some(StrandMmBins {
                    guide: guide_out.including_on_target,
                    passenger: pass_out.including_on_target,
                }),
                offtarget_only_genes: Some(StrandMmBins {
                    guide: guide_out.offtarget_only_genes,
                    passenger: pass_out.offtarget_only_genes,
                }),
                including_on_target_genes: Some(StrandMmBins {
                    guide: guide_out.including_on_target_genes,
                    passenger: pass_out.including_on_target_genes,
                }),
            }
        })
        .collect()
}

fn scan_transcript_t6b(
    store: &Store,
    tx: &TxMeta,
    tx_idx: u32,
    queries: &[PreparedQuery],
    filter: &TargetFilter,
    acc: &mut [(LocalStrand, LocalStrand)],
) {
    let n = tx.len as usize;
    if n < SEED_LEN {
        return;
    }
    let off = tx.offset as usize;
    let on_target = filter.is_on_target(tx);
    let gene_nonempty = !tx.gene_key.is_empty();

    for (qi, pq) in queries.iter().enumerate() {
        scan_pattern(
            store,
            tx,
            tx_idx,
            off,
            n,
            on_target,
            gene_nonempty,
            &pq.guide_pat,
            None,
            &mut acc[qi].0,
        );
        if let Some((_, pat)) = &pq.sense {
            scan_pattern(
                store,
                tx,
                tx_idx,
                off,
                n,
                on_target,
                gene_nonempty,
                pat,
                None,
                &mut acc[qi].1,
            );
        }
    }
}

fn scan_transcript_sidirect(
    store: &Store,
    tx: &TxMeta,
    tx_idx: u32,
    queries: &[PreparedQuery],
    filter: &TargetFilter,
    acc: &mut [(LocalSidirectStrand, LocalSidirectStrand)],
) {
    let n = tx.len as usize;
    if n < 19 {
        return;
    }
    let off = tx.offset as usize;
    let on_target = filter.is_on_target(tx);
    let gene_nonempty = !tx.gene_key.is_empty();

    for (qi, pq) in queries.iter().enumerate() {
        scan_pattern_sidirect(
            store,
            tx,
            tx_idx,
            off,
            n,
            on_target,
            gene_nonempty,
            &pq.guide_pat,
            "guide",
            &mut acc[qi].0,
        );
        if let Some((_, pat)) = &pq.sense {
            scan_pattern_sidirect(
                store,
                tx,
                tx_idx,
                off,
                n,
                on_target,
                gene_nonempty,
                pat,
                "passenger",
                &mut acc[qi].1,
            );
        }
    }
}

fn scan_pattern(
    store: &Store,
    tx: &TxMeta,
    tx_idx: u32,
    off: usize,
    n: usize,
    on_target: bool,
    gene_nonempty: bool,
    pat: &OligoPattern,
    strand: Option<&str>,
    out: &mut LocalStrand,
) {
    let l = pat.len;
    if n >= l {
        let last = n - l;
        for pos in 0..=last {
            if tx.has_n && window_has_n(&store.nmask, off + pos, l) {
                continue;
            }
            let window = extract_packed(&store.packed, off + pos, l);
            let mm = hamming_packed(window, pat.site, l);
            if mm <= MAX_MISMATCH as u32 {
                out.add_near(
                    Hit {
                        accession: tx.accession.clone(),
                        gene: tx.gene.clone(),
                        mismatches: mm as u8,
                        position: (pos + 1) as u32,
                        on_target,
                        site: decode_window(&store.packed, off + pos, l),
                        strand: strand.map(str::to_string),
                    },
                    on_target,
                    tx_idx,
                    tx.gene_id,
                );
            }
        }
    }

    if !pat.scan_seed || n < SEED_LEN {
        return;
    }
    let last7 = n - SEED_LEN;
    let mut saw_seed = false;
    for pos in 0..=last7 {
        if tx.has_n && window_has_n(&store.nmask, off + pos, SEED_LEN) {
            continue;
        }
        let window = extract_packed(&store.packed, off + pos, SEED_LEN);
        if window == pat.seed {
            saw_seed = true;
            break;
        }
    }
    if saw_seed {
        out.add_seed(on_target, tx_idx, tx.gene_id, gene_nonempty);
    }
}

fn scan_pattern_sidirect(
    store: &Store,
    tx: &TxMeta,
    tx_idx: u32,
    off: usize,
    n: usize,
    on_target: bool,
    gene_nonempty: bool,
    pat: &OligoPattern,
    strand: &str,
    out: &mut LocalSidirectStrand,
) {
    let l = pat.len;
    if n < l {
        return;
    }
    let last = n - l;
    let mut best: Option<u8> = None;
    for pos in 0..=last {
        if tx.has_n && window_has_n(&store.nmask, off + pos, l) {
            continue;
        }
        let window = extract_packed(&store.packed, off + pos, l);
        let mm = hamming_packed(window, pat.site, l);
        if mm <= MAX_MISMATCH as u32 {
            let mm_u8 = mm as u8;
            best = Some(best.map_or(mm_u8, |b| b.min(mm_u8)));
            push_hit(
                &mut out.hits,
                Hit {
                    accession: tx.accession.clone(),
                    gene: tx.gene.clone(),
                    mismatches: mm_u8,
                    position: (pos + 1) as u32,
                    on_target,
                    site: decode_window(&store.packed, off + pos, l),
                    strand: Some(strand.to_string()),
                },
            );
        }
    }
    if let Some(min_mm) = best {
        out.add_tx(min_mm, on_target, tx_idx, tx.gene_id, gene_nonempty);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::store::Store;
    use std::io::Write;

    fn fixture_store() -> Store {
        // guide ATAAACTCCAGGCCTATGAGG
        // site  CCTCATAGGCCTGGAGTTTAT
        // seed  TAAACTC → rc GAGTTTA
        let mut fa = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fa,
            ">NM_000001.1 Homo sapiens dummy target (DUMMY), mRNA\n\
             GGCCTCATAGGCCTGGAGTTTATGG\n\
             >NM_000002.1 Homo sapiens one mismatch (OFF1), mRNA\n\
             GGCCTCATAGGCCTGGAGTTTAAGG\n\
             >NM_000003.1 Homo sapiens two mismatch (OFF2), mRNA\n\
             AAAAAAAAAACCTCATAGGCCTGGAGTTTGG\n\
             >NM_000006.1 Homo sapiens three mismatch (OFF3), mRNA\n\
             AAAAAAAAAACCTCATAGGCCTGGAGTTGGG\n\
             >NR_000004.1 Homo sapiens seed only (SEEDG), long non-coding RNA\n\
             TTTTGAGTTTACCCCCCCCCCCCCC\n\
             >NM_000005.1 Homo sapiens perfect off (OFFP), mRNA\n\
             AAAACCTCATAGGCCTGGAGTTTATCCC"
        )
        .unwrap();
        fa.flush().unwrap();
        Store::from_fastas(&[fa.path().to_path_buf()], "test".into()).unwrap()
    }

    const GUIDE: &str = "ATAAACTCCAGGCCTATGAGG";

    #[test]
    fn planted_hits_and_on_target() {
        let store = fixture_store();
        let items = [SearchItem {
            id: "q1".into(),
            guide: GUIDE.into(),
            sense: None,
        }];
        let filter = TargetFilter::new(Some("DUMMY"), &[]);
        let r = search_batch(&store, &items, &filter, ScanMode::T6b { scan_sense: false });
        let r = &r[0];
        assert_eq!(r.on_target.transcripts, 1);
        assert_eq!(r.on_target.genes, 1);
        let off = r.offtarget.as_ref().unwrap();
        assert!(off.perfect >= 1, "{:?}", off);
        assert!(off.mismatch_1 >= 1, "{:?}", off);
        assert!(off.mismatch_2 >= 1, "{:?}", off);
        assert!(off.mismatch_3 >= 1, "{:?}", off);
        assert!(off.seed_transcripts >= 1);
        assert_eq!(r.specificity_label.as_deref(), Some("低"));
        assert!(r.hits.iter().any(|h| h.on_target && h.gene == "DUMMY"));
        assert!(r.hits.iter().any(|h| !h.on_target && h.gene == "OFFP"));
    }

    #[test]
    fn no_filter_counts_dummy_as_offtarget() {
        let store = fixture_store();
        let items = [SearchItem {
            id: "q1".into(),
            guide: GUIDE.into(),
            sense: None,
        }];
        let filter = TargetFilter::new(None, &[]);
        let r = &search_batch(&store, &items, &filter, ScanMode::T6b { scan_sense: false })[0];
        assert_eq!(r.on_target.transcripts, 0);
        assert!(r.offtarget.as_ref().unwrap().perfect >= 2);
        assert_eq!(r.specificity_label.as_deref(), Some("低"));
    }

    #[test]
    fn high_label_random() {
        let store = fixture_store();
        let items = [SearchItem {
            id: "rand".into(),
            guide: "GGGGGGGGGGGGGGGGGGG".into(),
            sense: None,
        }];
        let r = &search_batch(
            &store,
            &items,
            &TargetFilter::new(None, &[]),
            ScanMode::T6b { scan_sense: false },
        )[0];
        let off = r.offtarget.as_ref().unwrap();
        assert_eq!(off.perfect, 0);
        assert_eq!(off.mismatch_1, 0);
        if off.mismatch_2 == 0 && off.mismatch_3 == 0 && off.seed_transcripts == 0 {
            assert_eq!(r.specificity_label.as_deref(), Some("高"));
        }
    }

    #[test]
    fn seed_only_is_medium() {
        let store = fixture_store();
        let items = [SearchItem {
            id: "q1".into(),
            guide: GUIDE.into(),
            sense: None,
        }];
        let accs = vec![
            "NM_000001.1".into(),
            "NM_000002.1".into(),
            "NM_000003.1".into(),
            "NM_000005.1".into(),
            "NM_000006.1".into(),
        ];
        let filter = TargetFilter::new(None, &accs);
        let r = &search_batch(&store, &items, &filter, ScanMode::T6b { scan_sense: false })[0];
        let off = r.offtarget.as_ref().unwrap();
        assert_eq!(off.perfect, 0);
        assert_eq!(off.mismatch_1, 0);
        assert!(off.seed_transcripts >= 1);
        assert_eq!(r.specificity_label.as_deref(), Some("中"));
    }

    #[test]
    fn accession_filter() {
        let store = fixture_store();
        let items = [SearchItem {
            id: "q1".into(),
            guide: GUIDE.into(),
            sense: None,
        }];
        let filter = TargetFilter::new(None, &["NM_000001.1".into()]);
        let r = &search_batch(&store, &items, &filter, ScanMode::T6b { scan_sense: false })[0];
        assert_eq!(r.on_target.transcripts, 1);
        assert!(r
            .hits
            .iter()
            .any(|h| h.accession.starts_with("NM_000001") && h.on_target));
    }

    #[test]
    fn sidirect_both_strands_no_label() {
        let store = fixture_store();
        let items = [SearchItem {
            id: "q1".into(),
            guide: GUIDE.into(),
            sense: None,
        }];
        let filter = TargetFilter::new(Some("DUMMY"), &[]);
        let r = &search_batch(&store, &items, &filter, ScanMode::Sidirect)[0];
        assert_eq!(r.profile, "sidirect");
        assert!(r.specificity_label.is_none());
        assert!(r.offtarget.is_none());
        assert!(r.query_19mer.is_some());
        assert_eq!(
            r.query_19mer.as_ref().unwrap().guide,
            "TAAACTCCAGGCCTATGAG"
        );
        let mm = r.min_mismatch.as_ref().unwrap();
        // OFFP is perfect for full guide; 19-mer of guide may still hit
        assert!(mm.guide <= 3 || mm.guide == MIN_MISMATCH_NONE);
        assert!(r.passes_hide_less_specific.is_some());
        assert!(r.hits.iter().any(|h| h.strand.as_deref() == Some("guide")));
        // Passenger is always scanned; may have zero ≤3 hits in this tiny fixture.
        assert!(r.min_mismatch.as_ref().unwrap().passenger >= 2);
        assert_eq!(
            r.query_19mer.as_ref().unwrap().passenger.len(),
            19
        );
    }

    #[test]
    fn end_mismatch_affects_t6b_not_sidirect_19mer() {
        // Full 21-mer site with terminal mismatch vs 19-mer core that matches.
        // guide = ATAAACTCCAGGCCTATGAGG
        // 19mer =  TAAACTCCAGGCCTATGAG
        // site19 = CTCATAGGCCTGGAGTTT A  (rc of 19mer)
        // Plant a transcript that matches the 19-mer perfectly but differs at
        // the full-guide overhangs so full 21-mer Hamming = 2 (both ends).
        let mut fa = tempfile::NamedTempFile::new().unwrap();
        // Full guide site would be CCTCATAGGCCTGGAGTTTAT
        // 19-mer site (rc of TAAACTCCAGGCCTATGAG) = CTCATAGGCCTGGAGTTTA
        // Plant: flanking bases break full 21-mer but keep 19-mer perfect:
        // position: want window CTCATAGGCCTGGAGTTTA at offset making
        // a 21-nt with wrong ends: G + 19mer_site + G = GCTCATAGGCCTGGAGTTTAG
        // vs full site CCTCATAGGCCTGGAGTTTAT → mismatches at pos0 (G/C) and last (G/T) = 2
        writeln!(
            fa,
            ">NM_100001.1 Homo sapiens endmm (ENDMM), mRNA\n\
             AAAGCTCATAGGCCTGGAGTTTAGAAA\n\
             >NM_100002.1 Homo sapiens target (TARG), mRNA\n\
             AAACCTCATAGGCCTGGAGTTTATAAA"
        )
        .unwrap();
        fa.flush().unwrap();
        let store = Store::from_fastas(&[fa.path().to_path_buf()], "endmm".into()).unwrap();
        let items = [SearchItem {
            id: "q1".into(),
            guide: GUIDE.into(),
            sense: None,
        }];
        let filter = TargetFilter::new(Some("TARG"), &[]);

        let t6b = &search_batch(&store, &items, &filter, ScanMode::T6b { scan_sense: false })[0];
        let off = t6b.offtarget.as_ref().unwrap();
        assert!(
            off.mismatch_2 >= 1 || off.mismatch_1 >= 1 || off.perfect >= 1,
            "t6b should see ENDMM as near-match: {:?}",
            off
        );

        let sd = &search_batch(&store, &items, &filter, ScanMode::Sidirect)[0];
        // 19-mer perfect hit on ENDMM → offtarget_only guide[0] >= 1
        let bins = sd.offtarget_only.as_ref().unwrap();
        assert!(
            bins.guide[0] >= 1,
            "sidirect 19-mer should count ENDMM as perfect: {:?}",
            bins
        );
        assert_eq!(sd.min_mismatch.as_ref().unwrap().guide, 0);
    }

    #[test]
    fn sidirect_hide_when_min_ge_2() {
        // Guide with no near off-target in tiny empty-ish lib of polyA
        let mut fa = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            fa,
            ">NM_200001.1 Homo sapiens target (TARG), mRNA\n\
             {}\n\
             >NM_200002.1 Homo sapiens other (OTHER), mRNA\n\
             AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            // plant perfect on-target site for 19-mer of GUIDE
            // guide19 = TAAACTCCAGGCCTATGAG → site = CTCATAGGCCTGGAGTTTA
            format!("AAA{}AAA", "CTCATAGGCCTGGAGTTTA")
        )
        .unwrap();
        fa.flush().unwrap();
        let store = Store::from_fastas(&[fa.path().to_path_buf()], "hide".into()).unwrap();
        let items = [SearchItem {
            id: "q1".into(),
            guide: GUIDE.into(),
            sense: None,
        }];
        let r = &search_batch(
            &store,
            &items,
            &TargetFilter::new(Some("TARG"), &[]),
            ScanMode::Sidirect,
        )[0];
        let mm = r.min_mismatch.as_ref().unwrap();
        assert_eq!(mm.guide, MIN_MISMATCH_NONE);
        assert_eq!(mm.passenger, MIN_MISMATCH_NONE);
        assert_eq!(r.passes_hide_less_specific, Some(true));
    }
}
