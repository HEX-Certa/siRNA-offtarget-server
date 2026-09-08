use std::path::Path;

use redb::{Database, TableDefinition};

use crate::types::CachedResult;

const CACHE_TABLE: TableDefinition<&str, &[u8]> = TableDefinition::new("cache");
const META_TABLE: TableDefinition<&str, &str> = TableDefinition::new("meta");

const META_FP: &str = "fingerprint";
const META_SCHEMA: &str = "cache_schema";
/// Bump to invalidate all cached query results (no legacy key/value support).
pub const CACHE_SCHEMA: &str = "sidirect-api-1";

#[derive(Debug, thiserror::Error)]
pub enum CacheError {
    #[error("redb: {0}")]
    Redb(#[from] redb::Error),
    #[error("redb db: {0}")]
    Db(String),
}

impl From<redb::DatabaseError> for CacheError {
    fn from(e: redb::DatabaseError) -> Self {
        CacheError::Redb(e.into())
    }
}

impl From<redb::TransactionError> for CacheError {
    fn from(e: redb::TransactionError) -> Self {
        CacheError::Redb(e.into())
    }
}

impl From<redb::TableError> for CacheError {
    fn from(e: redb::TableError) -> Self {
        CacheError::Redb(e.into())
    }
}

impl From<redb::StorageError> for CacheError {
    fn from(e: redb::StorageError) -> Self {
        CacheError::Redb(e.into())
    }
}

impl From<redb::CommitError> for CacheError {
    fn from(e: redb::CommitError) -> Self {
        CacheError::Redb(e.into())
    }
}

pub struct Cache {
    db: Database,
}

impl Cache {
    pub fn open(path: &Path) -> Result<Self, CacheError> {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| CacheError::Db(e.to_string()))?;
        }
        let db = Database::create(path)?;
        let w = db.begin_write()?;
        {
            let _ = w.open_table(CACHE_TABLE)?;
            let _ = w.open_table(META_TABLE)?;
        }
        w.commit()?;
        Ok(Self { db })
    }

    /// Clear cache when FASTA fingerprint **or** API cache schema changes.
    pub fn ensure_fingerprint(&self, fp: &str) -> Result<(), CacheError> {
        let (stored_fp, stored_schema) = {
            let r = self.db.begin_read()?;
            let table = r.open_table(META_TABLE)?;
            let fp = table.get(META_FP)?.map(|v| v.value().to_string());
            let schema = table.get(META_SCHEMA)?.map(|v| v.value().to_string());
            (fp, schema)
        };
        if stored_fp.as_deref() == Some(fp) && stored_schema.as_deref() == Some(CACHE_SCHEMA) {
            return Ok(());
        }
        let w = self.db.begin_write()?;
        {
            w.delete_table(CACHE_TABLE)?;
            let _cache = w.open_table(CACHE_TABLE)?;
            let mut meta = w.open_table(META_TABLE)?;
            meta.insert(META_FP, fp)?;
            meta.insert(META_SCHEMA, CACHE_SCHEMA)?;
        }
        w.commit()?;
        Ok(())
    }

    pub fn get(&self, key: &str) -> Result<Option<CachedResult>, CacheError> {
        let r = self.db.begin_read()?;
        let table = r.open_table(CACHE_TABLE)?;
        match table.get(key)? {
            None => Ok(None),
            Some(v) => {
                let parsed =
                    serde_json::from_slice(v.value()).map_err(|e| CacheError::Db(e.to_string()))?;
                Ok(Some(parsed))
            }
        }
    }

    pub fn put(&self, key: &str, value: &CachedResult) -> Result<(), CacheError> {
        let bytes = serde_json::to_vec(value).map_err(|e| CacheError::Db(e.to_string()))?;
        let w = self.db.begin_write()?;
        {
            let mut table = w.open_table(CACHE_TABLE)?;
            table.insert(key, bytes.as_slice())?;
        }
        w.commit()?;
        Ok(())
    }
}

pub fn cache_key(
    lib_fp: &str,
    profile: &str,
    target_gene: &str,
    accessions_csv: &str,
    scan_sense: bool,
    guide: &str,
    sense: Option<&str>,
) -> String {
    let sense = sense.unwrap_or("");
    format!(
        "{lib_fp}|{profile}|{target_gene}|{accessions_csv}|{}|{guide}|{sense}",
        u8::from(scan_sense)
    )
}

pub fn normalize_target_gene(gene: Option<&str>) -> String {
    gene.map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| s.to_ascii_uppercase())
        .unwrap_or_default()
}

pub fn normalize_accessions_csv(accs: &[String]) -> String {
    let mut v: Vec<String> = accs
        .iter()
        .map(|s| s.trim().to_ascii_uppercase())
        .filter(|s| !s.is_empty())
        .collect();
    v.sort();
    v.dedup();
    v.join(",")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn key_stable_and_order_independent_accessions() {
        let a = cache_key(
            "fp",
            "t6b",
            "PCSK9",
            &normalize_accessions_csv(&["NM_2".into(), "NM_1".into()]),
            false,
            "ATAAACTCCAGGCCTATGAGG",
            None,
        );
        let b = cache_key(
            "fp",
            "t6b",
            "PCSK9",
            &normalize_accessions_csv(&["nm_1".into(), "NM_2".into()]),
            false,
            "ATAAACTCCAGGCCTATGAGG",
            None,
        );
        assert_eq!(a, b);
        assert!(a.starts_with("fp|t6b|PCSK9|NM_1,NM_2|0|"));
    }

    #[test]
    fn profiles_isolated_in_key() {
        let a = cache_key("fp", "t6b", "", "", false, "ATAAACTCCAGGCCTATGAGG", None);
        let b = cache_key(
            "fp",
            "sidirect",
            "",
            "",
            false,
            "ATAAACTCCAGGCCTATGAGG",
            None,
        );
        assert_ne!(a, b);
    }

    fn sample_cached() -> CachedResult {
        CachedResult {
            id: "x".into(),
            guide: "ATAAACTCCAGGCCTATGAGG".into(),
            length: 21,
            profile: "t6b".into(),
            specificity_label: Some("高".into()),
            on_target: Default::default(),
            offtarget: Some(Default::default()),
            hits: vec![],
            sense: None,
            query_19mer: None,
            min_mismatch: None,
            passes_hide_less_specific: None,
            offtarget_only: None,
            including_on_target: None,
            offtarget_only_genes: None,
            including_on_target_genes: None,
        }
    }

    #[test]
    fn put_get_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let cache = Cache::open(&dir.path().join("c.redb")).unwrap();
        cache.ensure_fingerprint("fp").unwrap();
        let val = sample_cached();
        let key = cache_key("fp", "t6b", "", "", false, &val.guide, None);
        cache.put(&key, &val).unwrap();
        let got = cache.get(&key).unwrap().unwrap();
        assert_eq!(got.id, "x");
        assert_eq!(got.specificity_label.as_deref(), Some("高"));
    }

    #[test]
    fn fingerprint_change_clears() {
        let dir = tempfile::tempdir().unwrap();
        let cache = Cache::open(&dir.path().join("c.redb")).unwrap();
        cache.ensure_fingerprint("fp1").unwrap();
        let val = sample_cached();
        let key = cache_key("fp1", "t6b", "", "", false, &val.guide, None);
        cache.put(&key, &val).unwrap();
        cache.ensure_fingerprint("fp2").unwrap();
        assert!(cache.get(&key).unwrap().is_none());
    }
}
