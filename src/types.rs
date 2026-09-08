use serde::{Deserialize, Serialize};

pub const DEFAULT_MAX_HITS: usize = 50;
pub const INTERNAL_HIT_CAP: usize = 200;
pub const DEFAULT_MAX_BATCH: usize = 50;

/// No non-target 0–3 hit ⇒ min mismatch reported as 4 (siDirect-compatible).
pub const MIN_MISMATCH_NONE: u8 = 4;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Profile {
    T6b,
    Sidirect,
}

impl Profile {
    pub fn as_str(self) -> &'static str {
        match self {
            Profile::T6b => "t6b",
            Profile::Sidirect => "sidirect",
        }
    }

    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "t6b" => Some(Profile::T6b),
            "sidirect" => Some(Profile::Sidirect),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct CheckRequest {
    #[serde(default)]
    pub target_gene: Option<String>,
    #[serde(default)]
    pub target_accessions: Option<Vec<String>>,
    pub queries: Vec<QueryIn>,
    #[serde(default)]
    pub scan_sense: bool,
    #[serde(default = "default_max_hits")]
    pub max_hits_per_query: usize,
    /// `"t6b"` (default) or `"sidirect"`. Parsed in the handler for a clear 400.
    #[serde(default = "default_profile_str")]
    pub profile: String,
}

fn default_max_hits() -> usize {
    DEFAULT_MAX_HITS
}

fn default_profile_str() -> String {
    "t6b".into()
}

#[derive(Debug, Clone, Deserialize)]
pub struct QueryIn {
    pub id: Option<String>,
    pub guide: String,
    pub sense: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CheckResponse {
    pub results: Vec<QueryResult>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryResult {
    pub id: String,
    pub guide: String,
    pub length: usize,
    pub cached: bool,
    pub profile: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub specificity_label: Option<String>,
    pub on_target: OnTargetCounts,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub offtarget: Option<OfftargetCounts>,
    pub hits: Vec<Hit>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sense: Option<StrandResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub query_19mer: Option<Query19mer>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_mismatch: Option<StrandMinMismatch>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub passes_hide_less_specific: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub offtarget_only: Option<StrandMmBins>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub including_on_target: Option<StrandMmBins>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub offtarget_only_genes: Option<StrandMmBins>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub including_on_target_genes: Option<StrandMmBins>,
}

/// Cached body: same as `QueryResult` minus `cached` (filled at serve time).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CachedResult {
    pub id: String,
    pub guide: String,
    pub length: usize,
    pub profile: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub specificity_label: Option<String>,
    pub on_target: OnTargetCounts,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub offtarget: Option<OfftargetCounts>,
    pub hits: Vec<Hit>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sense: Option<StrandResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub query_19mer: Option<Query19mer>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_mismatch: Option<StrandMinMismatch>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub passes_hide_less_specific: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub offtarget_only: Option<StrandMmBins>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub including_on_target: Option<StrandMmBins>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub offtarget_only_genes: Option<StrandMmBins>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub including_on_target_genes: Option<StrandMmBins>,
}

impl CachedResult {
    pub fn into_query_result(self, cached: bool, max_hits: usize) -> QueryResult {
        QueryResult {
            id: self.id,
            guide: self.guide,
            length: self.length,
            cached,
            profile: self.profile,
            specificity_label: self.specificity_label,
            on_target: self.on_target,
            offtarget: self.offtarget,
            hits: truncate_hits(self.hits, max_hits),
            sense: self.sense.map(|mut s| {
                s.hits = truncate_hits(s.hits, max_hits);
                s
            }),
            query_19mer: self.query_19mer,
            min_mismatch: self.min_mismatch,
            passes_hide_less_specific: self.passes_hide_less_specific,
            offtarget_only: self.offtarget_only,
            including_on_target: self.including_on_target,
            offtarget_only_genes: self.offtarget_only_genes,
            including_on_target_genes: self.including_on_target_genes,
        }
    }

    pub fn from_query_result(r: &QueryResult) -> Self {
        Self {
            id: r.id.clone(),
            guide: r.guide.clone(),
            length: r.length,
            profile: r.profile.clone(),
            specificity_label: r.specificity_label.clone(),
            on_target: r.on_target.clone(),
            offtarget: r.offtarget.clone(),
            hits: r.hits.clone(),
            sense: r.sense.clone(),
            query_19mer: r.query_19mer.clone(),
            min_mismatch: r.min_mismatch.clone(),
            passes_hide_less_specific: r.passes_hide_less_specific,
            offtarget_only: r.offtarget_only.clone(),
            including_on_target: r.including_on_target.clone(),
            offtarget_only_genes: r.offtarget_only_genes.clone(),
            including_on_target_genes: r.including_on_target_genes.clone(),
        }
    }
}

fn truncate_hits(mut hits: Vec<Hit>, max_hits: usize) -> Vec<Hit> {
    if hits.len() > max_hits {
        hits.truncate(max_hits);
    }
    hits
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StrandResult {
    pub sequence: String,
    pub specificity_label: String,
    pub on_target: OnTargetCounts,
    pub offtarget: OfftargetCounts,
    pub hits: Vec<Hit>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct OnTargetCounts {
    pub transcripts: u64,
    pub genes: u64,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct OfftargetCounts {
    pub perfect: u64,
    pub mismatch_1: u64,
    pub mismatch_2: u64,
    pub mismatch_3: u64,
    pub seed_transcripts: u64,
    pub seed_genes: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Query19mer {
    pub guide: String,
    pub passenger: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StrandMinMismatch {
    pub guide: u8,
    pub passenger: u8,
}

/// Per-strand mismatch bins `[n0, n1, n2, n3]` (transcript or gene units).
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct StrandMmBins {
    pub guide: [u64; 4],
    pub passenger: [u64; 4],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Hit {
    pub accession: String,
    pub gene: String,
    pub mismatches: u8,
    pub position: u32,
    pub on_target: bool,
    pub site: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct HealthResponse {
    pub status: &'static str,
}

#[derive(Debug, Clone, Serialize)]
pub struct ProfileInfo {
    pub name: &'static str,
    pub window: &'static str,
    pub strands: &'static str,
}

#[derive(Debug, Clone, Serialize)]
pub struct DbInfo {
    pub transcripts: u64,
    pub bases: u64,
    pub shards: Vec<String>,
    pub fingerprint: String,
    pub indexed: bool,
    pub includes_xm_xr: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub indexed_at: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub refseq_release: Option<String>,
    pub profiles: Vec<ProfileInfo>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ErrorBody {
    pub error: String,
}
