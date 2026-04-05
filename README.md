# NeoKCT
 
**Neo K-mer Count Table** a reference-free, memory-efficient framework for transcriptomic profiling at scale.
 
---
 
## Motivation
 
Reference-based analyses discard non-aligning reads, overlooking biologically important elements such as alternative open reading frames and aberrantly expressed tumor-specific antigens (aeTSAs). Systematic discovery of these sequences requires reference-free methods that scale to modern transcriptomic datasets comprising hundreds to thousands of RNA-Seq samples.
 
NeoKCT is built around a memory-efficient k-mer count table that holds billions of k-mers and their per-sample counts simultaneously, supports parallel traversal across large sample cohorts, and enables rapid reconstruction of transcriptomic profiles. By joining a table of sequenced aeTSA peptides with a table of RNA-Seq samples, one can find the origin of each peptide across thousands of samples entirely without reference alignment, enabling unbiased discovery of non-canonical transcripts for immunopeptidomics and transcriptomics applications.
 
---
 
## How It Works
 
```
RNA-Seq FASTQ files
        │
        ▼
  k-merize + count          (JelloFish.jl parallel, per-chunk)
        │  DNA k-mers → amino acid k-mers (5-bit encoding)
        ▼
  per-sample hash table
        │
        ▼
  push into NeoKCT           (NeoKCT.jl)
        │  CSR layout, bit-packed counts, prefix-indexed binary search
        ▼
  multi-sample KCT           ──► query / join / serialize
```
 
Each sample is processed into a k-mer count hash table, then merged into a single `NeoKCT` data structure that grows incrementally. The table can be periodically collapsed (deduplication of count words) and saved to disk in a versioned binary format.
 
---
 
## Key Technical Features
 
- **Bit-packed count storage** (`PackedArray`): variable-width values packed into fixed-size words, dramatically reducing memory footprint compared to standard arrays.
- **Compressed Sparse Row (CSR) layout**: k-mer sequences, offset arrays, and count word IDs are stored in three flat vectors, minimizing per-entry allocation overhead.
- **Prefix-indexed binary search**: a 4-symbol prefix index (20-bit for amino acids) partitions the sorted k-mer list into ~160,000 buckets, accelerating lookup.
- **Parallel k-merization**: multi-threaded chunk-based processing of FASTQ files with parallel hash-table merging (`jello_superthreaded_hash`).
- **Parallel merge-sort**: task-based divide-and-conquer sort used when integrating new samples into the table.
- **DNA → amino acid translation**: nucleotide k-mers are translated to amino acid k-mers using a custom 5-bit alphabet, enabling peptide-level queries.
- **Versioned binary serialization**: tables are written/read in a compact binary format (current: v1.2) with metadata headers.
 
---
 
## Requirements
 
- **Julia** 1.12+
 
Key dependencies (see `Project.toml`):
 
| Package | Purpose |
|---|---|
| `BioSequences` / `BioSymbols` | Biological sequence types |
| `Kmers` | K-mer representation and encoding |
| `GZip` | Compressed FASTQ support |
| `EzXML` | mzid (proteomics) file parsing |
| `CairoMakie` | Benchmarking plots |
| `ProgressMeter` | Progress bars |
 
Install dependencies from the repo root:
 
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
 
---
 
## Quick Start
 
```julia
include("NeoKCT.jl")
 
# Build a KCT from a single FASTQ file (K=30, amino-acid k-mers of length 10)
kct = build_kct("sample.fastq.gz")
 
# Build from multiple samples, saving checkpoints
samples = ["s1.fastq.gz", "s2.fastq.gz", "s3.fastq.gz"]
kct = build_kct(
    samples,
    30,                        # nucleotide k-mer length (must be divisible by 3)
    500_000;                   # chunk size for parallel processing
    word_size = UInt128,       # bit-packing word size
    collapse_every = 5,        # deduplicate count words every N samples
    save_at_samples = [10, 50, 100],
    save_path = "output/"
)
 
# Serialize / deserialize
write_kct(kct, "my_table.kct")
kct = read_kct("my_table.kct")
 
# Index a k-mer (returns Kmer => count_vector across samples)
i = findfirst(kct, kmer_bits)
kct[i]
```
 
---
 
## Project Layout
 
| File | Description |
|---|---|
| `NeoKCT.jl` | Core `NeoKCT` struct, `build_kct`, `push!`, binary search, sort |
| `PackedArray.jl` | Bit-packed variable-width array and deduplication |
| `JelloFish.jl` | Parallel k-merization and counting from FASTQ |
| `AAAlphabet.jl` | Custom 5-bit amino acid alphabet for `BioSequences` |
| `KCTLoader.jl` | Binary serialization / deserialization (v1.2) |
| `KCTBenchmarker.jl` | Memory and performance benchmarking with SVG plots |
| `BioParser.jl` | Unified reader for FASTQ, gzipped FASTQ, and mzid files |
| `parallel_sort.jl` | Task-based parallel merge-sort |
| `Project.toml` | Julia package manifest |
 
---
 
## Stability Notice
 
This project is under active development. The data structure layout, file format, and API are subject to significant change. No stable public API is guaranteed at this stage.

