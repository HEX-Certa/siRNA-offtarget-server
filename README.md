# siRNA Off-target Checker

对人 RefSeq RNA 做 siRNA 近完全匹配脱靶检查的内网服务。转录组建好后常驻内存，按批次查询，结果可缓存。

## 文档

- [接入指引](docs/接入指引.md) — HTTP 字段、示例、含靶对照表
- [原理与工作方式](docs/原理与工作方式.md) — 两种 profile、算法合同与文献出处

## 两种口径

| profile | 用途 |
|---------|------|
| `t6b`（默认） | 全长 guide Hamming 0–3 + seed（2–8）命中 + 高/中/低 |
| `sidirect` | oligo 2–20 的 19-mer、双侧、min mismatch / hide，便于和 siDirect 特异性检查对比 |

默认不是 siDirect。要比那一套数字，请显式传 `"profile": "sidirect"`。

## 快速试一下

```bash
BASE=http://127.0.0.1:8080   # 按部署修改

curl -sS "$BASE/health"
curl -sS "$BASE/v1/db/info"

curl -sS -X POST "$BASE/v1/offtarget/check" \
  -H 'Content-Type: application/json' \
  -d '{
    "profile": "sidirect",
    "target_gene": "PCSK9",
    "queries": [{"guide": "ATAAACTCCAGGCCTATGAGG"}]
  }'
```

## 本地跑

### 1. 下载 RefSeq RNA

需要 `data/human.1.rna.fna` … `human.16.rna.fna`（体积较大，已在 `.gitignore` 中忽略）。从 NCBI 拉取并解压：

```bash
./data/download-refseq-rna.sh
```

脚本默认下载分片 1–16，已存在的 `.fna` / `.gz` 会跳过。可覆盖范围，例如：

```bash
START=1 MAX=16 ./data/download-refseq-rna.sh
```

源地址默认是 `https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot`，可用环境变量 `BASE_URL` 覆盖。

### 2. 编译并启动

首次启动会解析 FASTA 并写出 `data/index/transcriptome.bin`（同样被 gitignore）。

```bash
cargo build --release
DATA_DIR=data LISTEN_ADDR=0.0.0.0:8080 ./target/release/sirna-offtarget-checker
```

常用环境变量：`DATA_DIR`、`INDEX_DIR`、`REDB_PATH`、`LISTEN_ADDR`、`MAX_BATCH`、`REFSEQ_RELEASE`。

```bash
cargo test
```

## 范围

只做全库近匹配（及 T6b 的 seed 计数）。不做候选设计、Ui-Tei 功效规则、seed Tm。无鉴权，按内网使用。
