use std::sync::Arc;

use axum::extract::State;
use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use axum::routing::{get, post};
use axum::{Json, Router};
use tower_http::trace::TraceLayer;

use crate::cache::{cache_key, normalize_accessions_csv, normalize_target_gene, Cache};
use crate::search::{search_batch, ScanMode, SearchItem, TargetFilter};
use crate::seq::{normalize, revcomp, validate_oligo, SeqError, MAX_OLIGO, MIN_OLIGO};
use crate::store::Store;
use crate::types::{
    CachedResult, CheckRequest, CheckResponse, DbInfo, ErrorBody, HealthResponse, Profile,
    ProfileInfo, QueryResult, INTERNAL_HIT_CAP,
};

#[derive(Clone)]
pub struct AppState {
    pub store: Arc<Store>,
    pub cache: Arc<Cache>,
    pub max_batch: usize,
    pub refseq_release: Option<String>,
}

pub fn router(state: AppState) -> Router {
    Router::new()
        .route("/health", get(health))
        .route("/v1/db/info", get(db_info))
        .route("/v1/offtarget/check", post(check))
        .layer(TraceLayer::new_for_http())
        .with_state(state)
}

async fn health() -> Json<HealthResponse> {
    Json(HealthResponse { status: "ok" })
}

fn profile_infos() -> Vec<ProfileInfo> {
    vec![
        ProfileInfo {
            name: "t6b",
            window: "full_guide",
            strands: "guide (optional sense via scan_sense)",
        },
        ProfileInfo {
            name: "sidirect",
            window: "oligo_2_20",
            strands: "both",
        },
    ]
}

async fn db_info(State(st): State<AppState>) -> Json<DbInfo> {
    Json(DbInfo {
        transcripts: st.store.transcripts(),
        bases: st.store.bases,
        shards: st.store.shards.clone(),
        fingerprint: st.store.fingerprint.clone(),
        indexed: true,
        includes_xm_xr: st.store.includes_xm_xr(),
        indexed_at: st.store.indexed_at,
        refseq_release: st.refseq_release.clone(),
        profiles: profile_infos(),
    })
}

struct ApiError {
    status: StatusCode,
    msg: String,
}

impl ApiError {
    fn bad(msg: impl Into<String>) -> Self {
        Self {
            status: StatusCode::BAD_REQUEST,
            msg: msg.into(),
        }
    }

    fn intern(msg: impl Into<String>) -> Self {
        Self {
            status: StatusCode::INTERNAL_SERVER_ERROR,
            msg: msg.into(),
        }
    }
}

impl IntoResponse for ApiError {
    fn into_response(self) -> Response {
        let body = Json(ErrorBody { error: self.msg });
        (self.status, body).into_response()
    }
}

async fn check(
    State(st): State<AppState>,
    Json(req): Json<CheckRequest>,
) -> Result<Json<CheckResponse>, ApiError> {
    if req.queries.is_empty() {
        return Err(ApiError::bad("queries must not be empty"));
    }
    if req.queries.len() > st.max_batch {
        return Err(ApiError::bad(format!(
            "at most {} siRNA sequences per request",
            st.max_batch
        )));
    }
    let max_hits = req.max_hits_per_query;
    if max_hits == 0 || max_hits > INTERNAL_HIT_CAP {
        return Err(ApiError::bad(format!(
            "max_hits_per_query must be 1–{INTERNAL_HIT_CAP}"
        )));
    }

    let profile = Profile::parse(&req.profile).ok_or_else(|| {
        ApiError::bad(format!(
            "profile must be \"t6b\" or \"sidirect\", got {:?}",
            req.profile
        ))
    })?;
    let mode = match profile {
        Profile::T6b => ScanMode::T6b {
            scan_sense: req.scan_sense,
        },
        Profile::Sidirect => ScanMode::Sidirect,
    };
    // sidirect ignores scan_sense for key stability (always both strands).
    let scan_sense_key = match profile {
        Profile::T6b => req.scan_sense,
        Profile::Sidirect => false,
    };

    let target_gene = normalize_target_gene(req.target_gene.as_deref());
    let accs = req.target_accessions.clone().unwrap_or_default();
    let acc_csv = normalize_accessions_csv(&accs);
    let filter = TargetFilter::new(
        if target_gene.is_empty() {
            None
        } else {
            Some(target_gene.as_str())
        },
        &accs,
    );

    let mut prepared: Vec<(String, SearchItem)> = Vec::with_capacity(req.queries.len());
    for (i, q) in req.queries.iter().enumerate() {
        let guide = normalize(&q.guide);
        validate_oligo(&guide).map_err(|e| seq_err("guide", e))?;
        let sense = match &q.sense {
            Some(s) => {
                let n = normalize(s);
                validate_oligo(&n).map_err(|e| seq_err("sense", e))?;
                Some(n)
            }
            None if matches!(profile, Profile::Sidirect) || req.scan_sense => {
                Some(revcomp(&guide))
            }
            None => None,
        };
        let id = q
            .id
            .clone()
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| format!("q{}", i + 1));
        let key = cache_key(
            &st.store.fingerprint,
            profile.as_str(),
            &target_gene,
            &acc_csv,
            scan_sense_key,
            &guide,
            sense.as_deref(),
        );
        prepared.push((key, SearchItem { id, guide, sense }));
    }

    let mut results: Vec<Option<QueryResult>> = vec![None; prepared.len()];
    let mut miss_items = Vec::new();
    let mut miss_keys = Vec::new();
    let mut miss_idx = Vec::new();

    for (i, (key, item)) in prepared.iter().enumerate() {
        match st.cache.get(key) {
            Ok(Some(cached)) => {
                let mut r = cached.into_query_result(true, max_hits);
                r.id = item.id.clone();
                results[i] = Some(r);
            }
            Ok(None) => {
                if let Some(prev) = miss_keys.iter().position(|k| k == key) {
                    miss_idx.push((i, prev));
                } else {
                    miss_keys.push(key.clone());
                    miss_idx.push((i, miss_items.len()));
                    miss_items.push(item.clone());
                }
            }
            Err(e) => return Err(ApiError::intern(e.to_string())),
        }
    }

    if !miss_items.is_empty() {
        let store = Arc::clone(&st.store);
        let computed =
            tokio::task::spawn_blocking(move || search_batch(&store, &miss_items, &filter, mode))
                .await
                .map_err(|e| ApiError::intern(e.to_string()))?;

        for (orig_i, miss_j) in miss_idx {
            if results[orig_i].is_some() {
                continue;
            }
            let mut r = computed[miss_j].clone();
            r.id = prepared[orig_i].1.id.clone();
            let cached = CachedResult::from_query_result(&r);
            if let Err(e) = st.cache.put(&prepared[orig_i].0, &cached) {
                tracing::warn!(error = %e, "cache put failed");
            }
            results[orig_i] = Some(r.into_truncated(max_hits));
        }
    }

    Ok(Json(CheckResponse {
        results: results.into_iter().map(|r| r.expect("filled")).collect(),
    }))
}

fn seq_err(field: &str, e: SeqError) -> ApiError {
    match e {
        SeqError::BadLength(n) => ApiError::bad(format!(
            "{field} must be {MIN_OLIGO}–{MAX_OLIGO} nt, got {n}"
        )),
        SeqError::BadBase(c) => ApiError::bad(format!("{field} contains invalid base {c:?}")),
        SeqError::Empty => ApiError::bad(format!("{field} is empty")),
    }
}

trait TruncateHits {
    fn into_truncated(self, max_hits: usize) -> QueryResult;
}

impl TruncateHits for QueryResult {
    fn into_truncated(self, max_hits: usize) -> QueryResult {
        CachedResult::from_query_result(&self).into_query_result(false, max_hits)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::store::Store;
    use axum::body::Body;
    use axum::http::{Request, StatusCode};
    use http_body_util::BodyExt;
    use tower::ServiceExt;

    fn test_app() -> Router {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("tiny.fna");
        std::fs::write(
            &fa,
            ">NM_000001.1 Homo sapiens dummy (DUMMY), mRNA\n\
             GGCCTCATAGGCCTGGAGTTTATGG\n",
        )
        .unwrap();
        let store = Store::from_fastas(&[fa], "api-fp".into()).unwrap();
        let cache = Cache::open(&dir.path().join("keep").join("c.redb")).unwrap();
        cache.ensure_fingerprint("api-fp").unwrap();
        std::mem::forget(dir);
        router(AppState {
            store: Arc::new(store),
            cache: Arc::new(cache),
            max_batch: 50,
            refseq_release: None,
        })
    }

    async fn send(app: Router, req: Request<Body>) -> (StatusCode, serde_json::Value) {
        let resp = app.oneshot(req).await.unwrap();
        let status = resp.status();
        let bytes = resp.into_body().collect().await.unwrap().to_bytes();
        let json: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
        (status, json)
    }

    #[tokio::test]
    async fn health_ok() {
        let (st, json) = send(
            test_app(),
            Request::builder()
                .uri("/health")
                .body(Body::empty())
                .unwrap(),
        )
        .await;
        assert_eq!(st, StatusCode::OK);
        assert_eq!(json["status"], "ok");
    }

    #[tokio::test]
    async fn db_info_has_profiles() {
        let (st, json) = send(
            test_app(),
            Request::builder()
                .uri("/v1/db/info")
                .body(Body::empty())
                .unwrap(),
        )
        .await;
        assert_eq!(st, StatusCode::OK);
        assert!(json["profiles"].as_array().unwrap().len() >= 2);
        assert!(json["includes_xm_xr"].is_boolean());
    }

    #[tokio::test]
    async fn rejects_over_batch() {
        let guides: Vec<_> = (0..51)
            .map(|i| {
                serde_json::json!({
                    "id": format!("q{i}"),
                    "guide": "ATAAACTCCAGGCCTATGAGG"
                })
            })
            .collect();
        let (st, json) = send(
            test_app(),
            Request::builder()
                .method("POST")
                .uri("/v1/offtarget/check")
                .header("content-type", "application/json")
                .body(Body::from(
                    serde_json::json!({ "queries": guides }).to_string(),
                ))
                .unwrap(),
        )
        .await;
        assert_eq!(st, StatusCode::BAD_REQUEST);
        assert!(json["error"].as_str().unwrap().contains("at most 50"));
    }

    #[tokio::test]
    async fn rejects_bad_length() {
        let (st, json) = send(
            test_app(),
            Request::builder()
                .method("POST")
                .uri("/v1/offtarget/check")
                .header("content-type", "application/json")
                .body(Body::from(
                    serde_json::json!({
                        "queries": [{ "id": "x", "guide": "ATAA" }]
                    })
                    .to_string(),
                ))
                .unwrap(),
        )
        .await;
        assert_eq!(st, StatusCode::BAD_REQUEST);
        assert!(json["error"].as_str().unwrap().contains("19"));
    }

    #[tokio::test]
    async fn rejects_bad_base() {
        let (st, json) = send(
            test_app(),
            Request::builder()
                .method("POST")
                .uri("/v1/offtarget/check")
                .header("content-type", "application/json")
                .body(Body::from(
                    serde_json::json!({
                        "queries": [{ "id": "x", "guide": "ATAAACTCCAGGCCTATGANN" }]
                    })
                    .to_string(),
                ))
                .unwrap(),
        )
        .await;
        assert_eq!(st, StatusCode::BAD_REQUEST);
        assert!(json["error"].as_str().unwrap().contains("invalid base"));
    }

    #[tokio::test]
    async fn rejects_bad_profile() {
        let (st, json) = send(
            test_app(),
            Request::builder()
                .method("POST")
                .uri("/v1/offtarget/check")
                .header("content-type", "application/json")
                .body(Body::from(
                    serde_json::json!({
                        "profile": "nope",
                        "queries": [{ "id": "x", "guide": "ATAAACTCCAGGCCTATGAGG" }]
                    })
                    .to_string(),
                ))
                .unwrap(),
        )
        .await;
        assert_eq!(st, StatusCode::BAD_REQUEST);
        assert!(json["error"].as_str().unwrap().contains("profile"));
    }

    #[tokio::test]
    async fn check_and_cache_t6b() {
        let app = test_app();
        let body = serde_json::json!({
            "target_gene": "DUMMY",
            "queries": [{
                "id": "PCSK9_21mer",
                "guide": "ATAAACTCCAGGCCTATGAGG"
            }]
        })
        .to_string();
        let (st, json) = send(
            app.clone(),
            Request::builder()
                .method("POST")
                .uri("/v1/offtarget/check")
                .header("content-type", "application/json")
                .body(Body::from(body.clone()))
                .unwrap(),
        )
        .await;
        assert_eq!(st, StatusCode::OK);
        assert_eq!(json["results"][0]["cached"], false);
        assert_eq!(json["results"][0]["profile"], "t6b");
        assert_eq!(json["results"][0]["on_target"]["transcripts"], 1);
        assert_eq!(json["results"][0]["id"], "PCSK9_21mer");
        assert!(json["results"][0]["specificity_label"].is_string());

        let (st2, json2) = send(
            app,
            Request::builder()
                .method("POST")
                .uri("/v1/offtarget/check")
                .header("content-type", "application/json")
                .body(Body::from(body))
                .unwrap(),
        )
        .await;
        assert_eq!(st2, StatusCode::OK);
        assert_eq!(json2["results"][0]["cached"], true);
        assert_eq!(
            json2["results"][0]["specificity_label"],
            json["results"][0]["specificity_label"]
        );
    }

    #[tokio::test]
    async fn sidirect_no_label_and_cache_isolated() {
        let app = test_app();
        let t6b_body = serde_json::json!({
            "target_gene": "DUMMY",
            "profile": "t6b",
            "queries": [{ "id": "q", "guide": "ATAAACTCCAGGCCTATGAGG" }]
        })
        .to_string();
        let sd_body = serde_json::json!({
            "target_gene": "DUMMY",
            "profile": "sidirect",
            "queries": [{ "id": "q", "guide": "ATAAACTCCAGGCCTATGAGG" }]
        })
        .to_string();

        let (_, t6b) = send(
            app.clone(),
            Request::builder()
                .method("POST")
                .uri("/v1/offtarget/check")
                .header("content-type", "application/json")
                .body(Body::from(t6b_body))
                .unwrap(),
        )
        .await;
        assert_eq!(t6b["results"][0]["profile"], "t6b");

        let (st, sd) = send(
            app.clone(),
            Request::builder()
                .method("POST")
                .uri("/v1/offtarget/check")
                .header("content-type", "application/json")
                .body(Body::from(sd_body.clone()))
                .unwrap(),
        )
        .await;
        assert_eq!(st, StatusCode::OK);
        assert_eq!(sd["results"][0]["profile"], "sidirect");
        assert!(sd["results"][0]["specificity_label"].is_null());
        assert!(sd["results"][0]["min_mismatch"].is_object());
        assert!(sd["results"][0]["passes_hide_less_specific"].is_boolean());
        assert_eq!(sd["results"][0]["cached"], false);

        let (_, sd2) = send(
            app,
            Request::builder()
                .method("POST")
                .uri("/v1/offtarget/check")
                .header("content-type", "application/json")
                .body(Body::from(sd_body))
                .unwrap(),
        )
        .await;
        assert_eq!(sd2["results"][0]["cached"], true);
    }
}
