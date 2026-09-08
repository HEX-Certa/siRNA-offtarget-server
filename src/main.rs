use std::env;
use std::fs;
use std::net::SocketAddr;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result};
use sirna_offtarget_checker::api::{router, AppState};
use sirna_offtarget_checker::cache::Cache;
use sirna_offtarget_checker::store::Store;
use sirna_offtarget_checker::types::DEFAULT_MAX_BATCH;
use tracing::info;

fn load_refseq_release(data_dir: &std::path::Path) -> Option<String> {
    if let Ok(v) = env::var("REFSEQ_RELEASE") {
        let t = v.trim().to_string();
        if !t.is_empty() {
            return Some(t);
        }
    }
    let path = data_dir.join("REFSEQ_RELEASE");
    fs::read_to_string(path).ok().and_then(|s| {
        let t = s.trim().to_string();
        if t.is_empty() {
            None
        } else {
            Some(t)
        }
    })
}

#[tokio::main]
async fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let data_dir = PathBuf::from(env::var("DATA_DIR").unwrap_or_else(|_| "data".into()));
    let index_dir = env::var("INDEX_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| data_dir.join("index"));
    let redb_path = env::var("REDB_PATH")
        .map(PathBuf::from)
        .unwrap_or_else(|_| index_dir.join("cache.redb"));
    let listen = env::var("LISTEN_ADDR").unwrap_or_else(|_| "0.0.0.0:8080".into());
    let max_batch = env::var("MAX_BATCH")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_MAX_BATCH);
    let refseq_release = load_refseq_release(&data_dir);

    info!(
        data_dir = %data_dir.display(),
        index_dir = %index_dir.display(),
        redb = %redb_path.display(),
        listen = %listen,
        "starting siRNA off-target service"
    );

    let store = Store::load_or_build(&data_dir, &index_dir).context("load transcriptome")?;
    let cache = Cache::open(&redb_path).context("open redb")?;
    cache
        .ensure_fingerprint(&store.fingerprint)
        .context("align cache fingerprint/schema")?;

    info!(
        transcripts = store.transcripts(),
        bases = store.bases,
        fingerprint = %store.fingerprint,
        includes_xm_xr = store.includes_xm_xr(),
        refseq_release = ?refseq_release,
        "transcriptome ready"
    );

    let app = router(AppState {
        store: Arc::new(store),
        cache: Arc::new(cache),
        max_batch,
        refseq_release,
    });

    let addr: SocketAddr = listen.parse().context("LISTEN_ADDR")?;
    let listener = tokio::net::TcpListener::bind(addr).await?;
    info!(%addr, "listening");
    axum::serve(listener, app).await?;
    Ok(())
}
