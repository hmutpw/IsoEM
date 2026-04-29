# IsoEM <img src="https://img.shields.io/badge/version-0.3.0-blue" align="right"/> <img src="https://img.shields.io/badge/R-%3E%3D4.1-brightgreen" align="right"/>

> EM-based transcript quantification for IsoQuant long-read RNA-seq output.  
> Supports bulk (single/multi-sample) and single-cell/spatial modes.  
> Designed for large files with memory-efficient chunk processing.

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
# From GitHub
devtools::install_github("hmutpw/IsoEM")

# From local tarball
install.packages("IsoEM_0.3.0.tar.gz", repos = NULL, type = "source")
```

**Dependencies:** `data.table`, `Matrix` (both on CRAN, auto-installed)

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
```

### Bulk — multi-sample (joint IsoQuant run)

```r
result <- run_isoem(
  counts_file    = "joint.transcript_model_reads.tsv.gz",
  gtf_file       = "joint.transcript_models.gtf",
  mode           = "bulk_multi",
  anno_file      = "read_to_sample.tsv",   # read_id | sample_id
  outdir         = "results/joint/",
  n_cores        = 8
)
# IsoEMDataset with count matrix (transcripts × samples)
```

### Single-cell — anno_file mode

```r
# Step 1: extract barcode/UMI from BAM (see below)
# Step 2: quantify
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
```

### Single-cell — regex mode (barcode/UMI embedded in read ID)

```r
# When read IDs contain barcode/UMI info (e.g. from BLAZE demultiplexing):
# READ_00001_bc=ACGTACGT_sq=TTTTAAAA
result <- run_sc_isoem(
  counts_file = "sample.transcript_model_reads.tsv.gz",
  gtf_file    = "sample.transcript_models.gtf",
  pattern     = ".*_bc=(?P<barcode>[ACGT]+)_sq=(?P<umi>[ACGT]+)",
  unit        = "umi",
  outdir      = "results/sc_regex/"
)
```

---

## Modular Step-by-Step Workflow

The modular API allows saving and reusing intermediate objects — useful when
tuning EM parameters on large datasets without re-processing raw files.

```r
# Step 1: validate inputs (no data loaded into memory)
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

# Step 3: run EM — rerun with different params without rebuilding EC
result <- run_em(ec, gtf_file = "sample.transcript_models.gtf",
                 max_iter = 500, tol = 1e-6, n_cores = 16)

# Step 4: write results
write_sc_isoem(result, outdir = "results/sc_sample/")
```

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
    -f valid_barcodes.txt   # optional whitelist

# Custom tags (e.g. STARsolo: CB:Z / UB:Z)
bash extract_barcodes_from_bam.sh \
    -i aligned.bam -o read_bc_umi.tsv.gz \
    -b CB:Z -u UB:Z
```

---

## Output Files

### Bulk (`write_isoem`)

```
outdir/
├── counts.tsv[.gz]      # transcript-level EM counts
└── qc.tsv[.gz]          # run QC summary
```

**`counts.tsv` columns:**

| Column | Description |
|--------|-------------|
| `transcript_id` | Transcript identifier |
| `sample_id` | Sample name |
| `em_count` | EM-estimated read count |
| `tpm` | Transcripts Per Million |
| `unique_count` | Uniquely assigned reads |
| `total_count` | All compatible reads |
| `certainty` | `unique / total` — quantification confidence (0–1) |
| `multimapping_rate` | `1 − certainty` |
| `gene_id` | Gene (from GTF) |
| `is_novel` | Novel IsoQuant transcript |

### Single-cell (`write_sc_isoem`)

```
outdir/
├── matrix/
│   ├── matrix.mtx[.gz]       # sparse count matrix (transcript × cell)
│   ├── features.tsv[.gz]     # transcript_id, gene_id, is_novel
│   └── barcodes.tsv[.gz]     # one barcode per line
└── cell_qc.tsv[.gz]          # per-cell QC metrics
```

Compatible with `Seurat::Read10X()` and `SingleCellExperiment`:

```r
# Seurat
library(Seurat)
seu <- Read10X("outdir/matrix/")
seu <- CreateSeuratObject(seu)

# SingleCellExperiment
library(BUSpaRse)
sce <- read10xCounts("outdir/matrix/")
```

**`cell_qc.tsv` columns:**

| Column | Description |
|--------|-------------|
| `barcode` | Cell barcode |
| `n_transcripts_detected` | Transcripts with EM count > 0.01 |
| `em_converged` | Whether EM converged |
| `em_n_iter` | Actual iterations |
| `n_ec` | Number of ECs |
| `n_unique_ec` | Single-transcript ECs |

---

## Toy Data Examples

All examples work out of the box with built-in toy data:

```r
library(IsoEM)

counts_f <- system.file("extdata", "toy_counts.tsv",       package = "IsoEM")
gtf_f    <- system.file("extdata", "toy.gtf",              package = "IsoEM")
bc_f     <- system.file("extdata", "toy_bc_umi.tsv",       package = "IsoEM")
multi_f  <- system.file("extdata", "toy_counts_multi.tsv", package = "IsoEM")
sid_f    <- system.file("extdata", "toy_sample_ids.tsv",   package = "IsoEM")
regex_f  <- system.file("extdata", "toy_counts_regex.tsv", package = "IsoEM")

# Bulk single-sample
run_isoem(counts_f, gtf_f, mode = "bulk_single",
          sample_id = "toy", outdir = tempdir())

# Bulk multi-sample
run_isoem(multi_f, gtf_f, mode = "bulk_multi",
          anno_file = sid_f, outdir = tempdir())

# Single-cell (anno_file)
run_sc_isoem(counts_f, gtf_f, anno_file = bc_f,
             unit = "umi", outdir = tempdir())

# Single-cell (regex)
run_sc_isoem(regex_f, gtf_f,
             pattern = ".*_bc=(?P<barcode>[A-Z0-9]+)_sq=(?P<umi>[ACGT]+)",
             unit = "umi", outdir = tempdir())
```

---

## Notes for TE Analysis

IsoEM was designed with **transposable element transcript quantification** in mind:

- **High multimapping:** TE loci (e.g. L1 LINE) show `certainty < 0.1` because
  a single read may be compatible with hundreds of near-identical copies.
  This is **expected** — EM distributes counts probabilistically.
- **EC compression:** Large ECs (100+ transcripts) from TE loci are handled
  efficiently; the EM algorithm scales with EC count, not read count.
- **Novel transcripts:** IsoQuant-assembled TE-derived isoforms are flagged as
  `is_novel = TRUE` in all outputs.
- **UMI mode for SC:** Always use `unit = "umi"` for single-cell data.
  EC construction uses `(barcode, UMI)` as the observational unit, correctly
  handling UMI collisions across barcodes.

---

## Memory Efficiency

IsoEM handles very large files through:

| Strategy | Detail |
|----------|--------|
| Lazy loading | `prepare_isoem()` validates paths without reading data |
| Integer mapping | Reads and transcripts converted to integers before EC construction |
| Chunk processing | `anno_file` / `counts_file` read in configurable chunks (`chunk_size`) |
| Early filtering | `anno_file` records not in `counts_file` discarded immediately |
| EC compression | EM operates on EC table (small) not raw reads (large) |
| `keep_temp` | Intermediate merged files deleted after EC construction by default |

For a typical spatial transcriptomics sample (45M reads, 120M barcode records):
estimated peak memory ~4 GB with default settings.

---

## Citation

If you use IsoEM in your research, please cite:

> hmutpw. *IsoEM: EM-based isoform quantification for IsoQuant long-read output.*
> GitHub: https://github.com/hmutpw/IsoEM (2025)

---

## License

MIT © [hmutpw](https://github.com/hmutpw/IsoEM)
