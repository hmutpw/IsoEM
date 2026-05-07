# IsoEM <img src="https://img.shields.io/badge/version-0.3.1-blue" align="right"/> <img src="https://img.shields.io/badge/R-%3E%3D4.1-brightgreen" align="right"/> <img src="https://img.shields.io/badge/license-MIT-lightgrey" align="right"/>

> EM-based transcript quantification for [IsoQuant](https://github.com/ablab/isoquant) long-read RNA-seq output.  
> Supports bulk (single/multi-sample) and single-cell/spatial modes.  
> Designed for large FFPE and single-cell datasets with memory-efficient chunk processing.

---

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
  - [Bulk — single sample](#bulk--single-sample)
  - [Bulk — multi-sample](#bulk--multi-sample)
  - [Single-cell — anno\_file mode](#single-cell--anno_file-mode)
  - [Single-cell — regex mode](#single-cell--regex-mode)
- [Modular Step-by-Step Workflow](#modular-step-by-step-workflow)
- [Output Files](#output-files)
  - [Bulk outputs](#bulk-write_isoem)
  - [Single-cell outputs](#single-cell-write_sc_isoem)
  - [EC table](#ec-table-ec_tabletsv)
  - [Sharing table](#sharing-table-sharing_tabletsv)
- [Extracting Barcode/UMI from BAM](#extracting-barcodesumi-from-bam)
- [Loading Results Downstream](#loading-results-downstream)
- [Notes for TE Analysis](#notes-for-te-analysis)
- [Memory Efficiency](#memory-efficiency)
- [Changelog](#changelog)
- [Citation](#citation)

---

## Overview

IsoEM post-processes [IsoQuant](https://github.com/ablab/isoquant) output using
**equivalence-class (EC) compressed EM** to produce statistically principled
transcript-level counts. The pipeline is fully modular:

```
prepare_isoem()   →  IsoEMInput    validate paths, no data loaded
      ↓
build_ec()        →  IsoEMEC       bulk: build ECs
build_sc_ec()     →  IsoEMEC       SC:   build ECs + UMI deduplication
      ↓
run_em()          →  IsoEMResult / IsoEMDataset / IsoEMSCResult
      ↓
write_isoem()     /  write_sc_isoem()    write results to disk
```

One-step wrappers `run_isoem()` and `run_sc_isoem()` chain all steps automatically.

---

## Installation

```r
# From GitHub (recommended)
if (!require("remotes")) install.packages("remotes")
remotes::install_github("hmutpw/IsoEM", upgrade = "never")

# From local tarball
install.packages("IsoEM_0.3.1.tar.gz", repos = NULL, type = "source")
```

**Dependencies:** `data.table`, `Matrix` (both on CRAN, installed automatically)

---

## Quick Start

### Bulk — single sample

```r
library(IsoEM)

result <- run_isoem(
  counts_file = "sample.transcript_model_reads.tsv.gz",
  gtf_file    = "sample.transcript_models.gtf",
  mode        = "bulk_single",
  sample_id   = "my_sample",
  outdir      = "results/my_sample/"
)

print(result)
# IsoEMResult
#   Sample             : my_sample
#   Transcripts        : 15234
#   Novel transcripts  : 412
#   Unique assign rate : 68.3%
#   Converged          : TRUE (87 iterations)

# Write results to disk
write_isoem(result, outdir = "results/my_sample/")
# results/my_sample/
# ├── counts.tsv
# ├── qc.tsv
# ├── ec_table.tsv       ← NEW in v0.3.1
# └── sharing_table.tsv  ← NEW in v0.3.1
```

### Bulk — multi-sample

```r
result <- run_isoem(
  counts_file = "joint.transcript_model_reads.tsv.gz",
  gtf_file    = "joint.transcript_models.gtf",
  mode        = "bulk_multi",
  anno_file   = "read_to_sample.tsv",   # read_id | sample_id
  outdir      = "results/joint/",
  n_cores     = 8
)
# IsoEMDataset with count matrix (transcripts × samples)

write_isoem(result, outdir = "results/joint/")
```

### Single-cell — anno_file mode

```r
result <- run_sc_isoem(
  counts_file  = "sample.transcript_model_reads.tsv.gz",
  gtf_file     = "sample.transcript_models.gtf",
  anno_file    = "read_bc_umi.tsv.gz",     # read_id | barcode | umi
  barcodes_use = "valid_barcodes.txt",     # optional whitelist
  unit         = "umi",
  outdir       = "results/sc_sample/",
  n_cores      = 16
)

print(result)
# IsoEMSCResult
#   Cells              : 5432
#   Transcripts        : 15234
#   Median UMI / cell  : 312
#   Cells converged    : 5421 / 5432

write_sc_isoem(result, outdir = "results/sc_sample/")
# results/sc_sample/
# ├── matrix/
# │   ├── matrix.mtx
# │   ├── features.tsv
# │   └── barcodes.tsv
# ├── cell_qc.tsv
# ├── ec_table.tsv       ← NEW in v0.3.1
# └── sharing_table.tsv  ← NEW in v0.3.1
```

### Single-cell — regex mode

For read IDs that contain barcode/UMI information (e.g. from BLAZE demultiplexing):

```
READ_00001_bc=ACGTACGT_sq=TTTTAAAA
```

```r
result <- run_sc_isoem(
  counts_file = "sample.transcript_model_reads.tsv.gz",
  gtf_file    = "sample.transcript_models.gtf",
  pattern     = ".*_bc=(?P<barcode>[ACGT]+)_sq=(?P<umi>[ACGT]+)",
  unit        = "umi",
  outdir      = "results/sc_regex/"
)
```

Alternatively, supply `bc_pattern` and `umi_pattern` as separate regex strings
if barcode and UMI are in different fields.

---

## Modular Step-by-Step Workflow

The modular API allows saving and reusing intermediate objects — useful for
tuning EM parameters on large datasets without re-processing raw files.

```r
# Step 1: validate inputs (no data loaded into memory yet)
input <- prepare_isoem(
  counts_file  = "sample.transcript_model_reads.tsv.gz",
  gtf_file     = "sample.transcript_models.gtf",
  mode         = "sc",
  anno_file    = "read_bc_umi.tsv.gz",
  barcodes_use = "valid_barcodes.txt",
  unit         = "umi",
  outdir       = "results/sc_sample/"
)
print(input)
# IsoEMInput
#   Mode         : sc
#   counts_file  : sample.transcript_model_reads.tsv.gz (2.3 GB, ~45M rows)
#   anno_file    : read_bc_umi.tsv.gz (8.1 GB, ~120M rows)
#   barcodes_use : 6234 barcodes
#   unit         : umi
#   Validated    : TRUE
#   Est. memory  : ~4.0 GB peak

# Step 2: build EC (memory-efficient, chunk processing)
ec <- build_sc_ec(input, n_cores = 16)
saveRDS(ec, "results/sc_sample/ec.rds")   # save for reuse

print(ec)
# IsoEMEC
#   Mode : sc  |  Unit : umi
#   Groups (cells)  : 6234
#   Transcripts     : 15234
#   Total ECs       : 1823456
#   Median EC/cell  : 287

# Step 3: run EM — can be re-run with different parameters without rebuilding EC
result <- run_em(
  ec       = ec,
  gtf_file = "sample.transcript_models.gtf",
  max_iter = 500,
  tol      = 1e-6,
  n_cores  = 16
)

# Step 4: write results
write_sc_isoem(
  result,
  outdir              = "results/sc_sample/",
  compress            = TRUE,           # write .tsv.gz
  write_ec_table      = TRUE,
  write_sharing_table = TRUE,
  min_sharing         = 0.5,            # only report pairs sharing ≥50% of reads
  max_ec_size_sharing = 50L             # skip repeat-heavy ECs > 50 transcripts
)
```

### Reloading a saved EC object

```r
# Load saved EC and re-run EM with different parameters
ec     <- readRDS("results/sc_sample/ec.rds")
result <- run_em(ec, gtf_file = "sample.transcript_models.gtf",
                 max_iter = 1000, tol = 1e-8, n_cores = 32)
write_sc_isoem(result, outdir = "results/sc_sample_v2/")
```

---

## Output Files

### Bulk (`write_isoem`)

```
outdir/
├── counts.tsv[.gz]         # transcript-level EM counts
├── qc.tsv[.gz]             # run QC summary
├── ec_table.tsv[.gz]       # equivalence class table    [v0.3.1]
└── sharing_table.tsv[.gz]  # transcript sharing pairs   [v0.3.1]

# Multi-sample: per-sample subdirs + combined matrices
outdir/
├── sampleA/
│   ├── counts.tsv
│   ├── qc.tsv
│   ├── ec_table.tsv
│   └── sharing_table.tsv
├── sampleB/
│   └── ...
└── count_matrix.tsv        # transcripts × samples
```

**`counts.tsv` columns:**

| Column              | Description                                        |
| ------------------- | -------------------------------------------------- |
| `transcript_id`     | Transcript identifier                              |
| `sample_id`         | Sample name                                        |
| `em_count`          | EM-estimated read count                            |
| `tpm`               | Transcripts Per Million                            |
| `unique_count`      | Uniquely assigned reads                            |
| `total_count`       | All compatible reads                               |
| `certainty`         | `unique / total` — quantification confidence (0–1) |
| `multimapping_rate` | `1 − certainty`                                    |
| `gene_id`           | Gene (from GTF)                                    |
| `is_novel`          | Novel IsoQuant transcript                          |

**`qc.tsv` columns:**

| Column                   | Description                      |
| ------------------------ | -------------------------------- |
| `sample_id`              | Sample name                      |
| `total_reads`            | Total reads assigned to any EC   |
| `unique_assignment_rate` | Fraction of reads in unique ECs  |
| `n_transcripts_detected` | Transcripts with EM count > 0.01 |
| `n_unique_ecs`           | Single-transcript ECs            |
| `n_multi_ecs`            | Multi-transcript ECs             |
| `n_iter`                 | EM iterations used               |
| `converged`              | Whether EM converged             |

---

### Single-cell (`write_sc_isoem`)

```
outdir/
├── matrix/
│   ├── matrix.mtx[.gz]       # sparse count matrix (transcript × cell)
│   ├── features.tsv[.gz]     # transcript_id, gene_id, is_novel
│   └── barcodes.tsv[.gz]     # one barcode per line
├── cell_qc.tsv[.gz]          # per-cell QC metrics
├── ec_table.tsv[.gz]         # equivalence class table    [v0.3.1]
└── sharing_table.tsv[.gz]    # transcript sharing pairs   [v0.3.1]
```

**`cell_qc.tsv` columns:**

| Column                   | Description                                               |
| ------------------------ | --------------------------------------------------------- |
| `barcode`                | Cell barcode                                              |
| `n_transcripts_detected` | Transcripts with EM count > 0.01                          |
| `em_converged`           | Whether EM converged for this cell                        |
| `em_n_iter`              | Actual iterations used                                    |
| `n_ec`                   | Number of ECs for this cell                               |
| `n_unique_ec`            | Single-transcript ECs (ec_size == 1)                      |
| `total_reads`            | Total UMIs assigned to this cell                          |
| `unique_reads`           | UMIs from unique-mapping ECs                              |
| `unique_read_frac`       | `unique_reads / total_reads` — mapping specificity (0--1) |

---

### EC table (`ec_table.tsv`)

Records every equivalence class and its associated transcripts. Available for
both bulk and SC modes.

**Bulk columns:**

| Column        | Description                                            |
| ------------- | ------------------------------------------------------ |
| `sample_id`   | Sample name                                            |
| `ec_id`       | Integer EC identifier (unique within sample)           |
| `transcripts` | Pipe-separated transcript names (e.g. `TX1\|TX2\|TX3`) |
| `count`       | Number of reads/UMIs in this EC                        |
| `ec_size`     | Number of compatible transcripts                       |
| `ec_type`     | `unique` (ec_size = 1) or `multi` (ec_size > 1)        |

**SC columns:** same, with `group_id` (barcode) instead of `sample_id`.

> **Note:** ECs with `ec_type = "multi"` represent multi-mapping reads that are
> probabilistically distributed across transcripts by the EM algorithm.
> High `ec_size` values at TE loci are expected and handled correctly.

---

### Sharing table (`sharing_table.tsv`)

Identifies transcript pairs that share a large fraction of their reads through
multi-mapping ECs. Useful for assessing quantification reliability, especially
at repetitive TE loci.

**Sharing fraction** is defined as:

```
sharing_fraction = shared_reads / min(total_reads_tx1, total_reads_tx2)
```

**Columns:**

| Column                   | Description                                             |
| ------------------------ | ------------------------------------------------------- |
| `sample_id` / `group_id` | Sample name or cell barcode                             |
| `transcript_1`           | First transcript                                        |
| `transcript_2`           | Second transcript                                       |
| `shared_reads`           | Reads in ECs containing both transcripts                |
| `total_reads_tx1`        | Total reads compatible with transcript 1                |
| `total_reads_tx2`        | Total reads compatible with transcript 2                |
| `sharing_fraction`       | `shared / min(total_tx1, total_tx2)`                    |
| `recommendation`         | `consider_merging` (≥0.9) or `ambiguous_quantification` |

**Control parameters** in `write_isoem()` / `write_sc_isoem()`:

```r
write_isoem(result, outdir,
  write_ec_table      = TRUE,   # write ec_table.tsv
  write_sharing_table = TRUE,   # write sharing_table.tsv
  min_sharing         = 0.5,    # only report pairs with sharing ≥ 0.5
  max_ec_size_sharing = 50L     # skip ECs with > 50 transcripts (TE repeat clusters)
)
```

Setting `max_ec_size_sharing` prevents combinatorial explosion at highly
repetitive loci (e.g. L1 LINEs with 500+ near-identical copies).
These large ECs are simply excluded from pairwise sharing computation.

---

## Extracting Barcode/UMI from BAM

IsoEM includes a shell script for extracting `read_id | barcode | umi`
from minimap2-aligned BAM files:

```bash
# Default tags: BC:Z (barcode), U8:Z (UMI)
bash $(Rscript -e "cat(system.file('scripts',
    'extract_barcodes_from_bam.sh', package='IsoEM'))") \
    -i aligned.bam \
    -o read_bc_umi.tsv.gz \
    -t 8 \
    -f valid_barcodes.txt   # optional whitelist filter

# Custom BAM tags (e.g. STARsolo output: CB:Z / UB:Z)
bash extract_barcodes_from_bam.sh \
    -i aligned.bam -o read_bc_umi.tsv.gz \
    -b CB:Z -u UB:Z
```

**Output format** (`read_bc_umi.tsv.gz`):

```
READ_00001    ACGTACGT    TTTTAAAA
READ_00002    TTGCAAGT    GCGCGCGC
...
```

3 columns, no header: `read_id | barcode | umi`

---

## Loading Results Downstream

### Seurat (single-cell)

```r
library(Seurat)

# Load sparse count matrix
seu <- Read10X("results/sc_sample/matrix/")
seu <- CreateSeuratObject(seu, project = "IsoEM")

# Cell QC already computed — add to metadata
cell_qc <- read.table("results/sc_sample/cell_qc.tsv",
                       header = TRUE, sep = "\t", row.names = 1)
seu <- AddMetaData(seu, cell_qc)

# Filter by EM convergence
seu <- subset(seu, em_converged == TRUE)
```

### SingleCellExperiment

```r
library(BUSpaRse)
sce <- read10xCounts("results/sc_sample/matrix/")
```

### DESeq2 (bulk)

```r
library(data.table)
counts_dt <- fread("results/joint/count_matrix.tsv")
count_mat <- as.matrix(counts_dt[, -1], rownames = counts_dt$transcript_id)
count_mat <- round(count_mat)   # DESeq2 requires integers

library(DESeq2)
dds <- DESeqDataSetFromMatrix(count_mat, colData = sample_info, design = ~condition)
```

### Inspecting EC and sharing tables

```r
library(data.table)

# EC table: how many multi-mapping ECs per sample?
ec <- fread("results/my_sample/ec_table.tsv")
ec[, .N, by = ec_type]
#    ec_type      N
# 1:  unique  12043
# 2:   multi   3191

# Sharing table: transcripts that cannot be distinguished
sh <- fread("results/my_sample/sharing_table.tsv")
sh[recommendation == "consider_merging"][order(-sharing_fraction)]
#    sample_id  transcript_1  transcript_2  sharing_fraction  recommendation
# 1: my_sample  TX_L1_001     TX_L1_002     0.97              consider_merging
# 2: my_sample  TX_ERV_003    TX_ERV_004    0.93              consider_merging
```

---

## Toy Data Examples

All examples work with built-in toy data:

```r
library(IsoEM)

counts_f <- system.file("extdata", "toy_counts.tsv",       package = "IsoEM")
gtf_f    <- system.file("extdata", "toy.gtf",              package = "IsoEM")
bc_f     <- system.file("extdata", "toy_bc_umi.tsv",       package = "IsoEM")
multi_f  <- system.file("extdata", "toy_counts_multi.tsv", package = "IsoEM")
sid_f    <- system.file("extdata", "toy_sample_ids.tsv",   package = "IsoEM")
regex_f  <- system.file("extdata", "toy_counts_regex.tsv", package = "IsoEM")

# Bulk single-sample
result <- run_isoem(counts_f, gtf_f, mode = "bulk_single",
                    sample_id = "toy", outdir = tempdir())
write_isoem(result, outdir = file.path(tempdir(), "toy_bulk"))

# Bulk multi-sample
result <- run_isoem(multi_f, gtf_f, mode = "bulk_multi",
                    anno_file = sid_f, outdir = tempdir())
write_isoem(result, outdir = file.path(tempdir(), "toy_multi"))

# Single-cell (anno_file)
result <- run_sc_isoem(counts_f, gtf_f, anno_file = bc_f,
                       unit = "umi", outdir = tempdir())
write_sc_isoem(result, outdir = file.path(tempdir(), "toy_sc"))

# Single-cell (regex)
result <- run_sc_isoem(regex_f, gtf_f,
                       pattern = ".*_bc=(?P<barcode>[A-Z0-9]+)_sq=(?P<umi>[ACGT]+)",
                       unit = "umi", outdir = tempdir())
```

---

## Notes for TE Analysis

IsoEM was designed with **transposable element transcript quantification** in mind:

- **High multimapping:** TE loci (e.g. L1 LINEs) show `certainty < 0.1` because
  a single read may be compatible with hundreds of near-identical copies.
  This is **expected** — EM distributes counts probabilistically across all
  compatible transcripts. Do not filter on certainty alone for TE studies.

- **EC compression:** Large ECs (100+ transcripts) from TE loci are handled
  efficiently; the EM algorithm scales with EC count, not read count.

- **Novel transcripts:** IsoQuant-assembled TE-derived isoforms are flagged
  `is_novel = TRUE` in all outputs.

- **UMI mode for SC:** Always use `unit = "umi"` for single-cell data.
  EC construction uses `(barcode, UMI)` as the observational unit, correctly
  handling UMI collisions across barcodes.

- **Sharing table for TE QC:** The `sharing_table.tsv` is particularly valuable
  for TE work — pairs with `sharing_fraction ≥ 0.9` at TE loci typically
  indicate subfamilies that cannot be distinguished at read level. Consider
  aggregating these to the subfamily level for downstream analysis.
  Use `max_ec_size_sharing` (default 50) to prevent the pairwise expansion from
  becoming intractable at highly repetitive loci.

---

## Memory Efficiency

IsoEM handles very large files through several strategies:

| Strategy         | Detail                                                                                  |
| ---------------- | --------------------------------------------------------------------------------------- |
| Lazy loading     | `prepare_isoem()` validates paths without reading data                                  |
| Integer mapping  | Reads and transcripts integerised before EC construction                                |
| Chunk processing | `anno_file` / `counts_file` read in configurable chunks (`chunk_size`, default 5M rows) |
| Early filtering  | `anno_file` records not in `counts_file` discarded immediately                          |
| EC compression   | EM operates on the EC table (small), not raw reads (large)                              |
| Temp cleanup     | Intermediate files deleted after EC construction by default (`keep_temp = FALSE`)       |

**Typical memory usage** for a spatial transcriptomics sample
(45M reads, 120M barcode records): **~4 GB peak** with default settings.

For very large datasets, reduce `chunk_size` to lower peak memory at the cost
of slightly longer I/O time:

```r
ec <- build_sc_ec(input, chunk_size = 1e6L, n_cores = 16)
```

---

## Changelog

### v0.3.1

- **New:** `ec_table.tsv` output — all equivalence classes with transcript
  membership, count, size, and type (`unique`/`multi`)
- **New:** `sharing_table.tsv` output — pairwise transcript sharing analysis
  from multi-mapping ECs, with `sharing_fraction` and `recommendation` columns
- **New parameters** in `write_isoem()` and `write_sc_isoem()`:
  `write_ec_table`, `write_sharing_table`, `min_sharing`, `max_ec_size_sharing`
- EC data is now stored in `IsoEMResult` and `IsoEMSCResult` objects and
  carried through from `build_ec()` / `build_sc_ec()` automatically
- Fully backward-compatible — existing code produces identical output by default

### v0.3.0

- Initial public release

---

## Citation

If you use IsoEM in your research, please cite:

> hmutpw. *IsoEM: EM-based isoform quantification for IsoQuant long-read output.*  
> GitHub: https://github.com/hmutpw/IsoEM (2025)

---

## License

MIT © [hmutpw](https://github.com/hmutpw/IsoEM)
