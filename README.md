# Automated-Autocycler-Pipeline
Run the full AutoCycler end-to-end pipeline over all reads in a folder — hands-free. This script discovers FASTQ/BAM files, filters, subsamples, assembles with multiple assemblers, clusters, trims, resolves, combines, polishes with Racon, and finally assesses quality with CheckM2.

- Supports inputs: .fastq, .fq, .fastq.gz, .fq.gz, and .bam (BAM is auto-converted to FASTQ)
- Caches results per step so you can safely resume/re-run
- Organizes outputs by step and isolate
- Chooses appropriate mapping preset for minimap2 based on read type
- Ends with a CheckM2 report

Suitable for high-throughput of ONT and PacBio long reads.

---

## Contents
- What this script does
- Requirements
- Quick start
- Configuration
- Output structure
- Resumability and caching
- Notes and tips
- Troubleshooting
- Citations and acknowledgements

---

## What this script does

Pipeline overview:
1) Filter reads (filtlong)  
2) Estimate genome size (AutoCycler helper) with fallback  
3) Subsample reads (AutoCycler)  
4) Assemble with multiple assemblers (AutoCycler helpers)  
5) Compress assemblies (AutoCycler compress)  
6) Cluster contigs (AutoCycler cluster)  
7) Trim and resolve clusters (AutoCycler trim/resolve)  
8) Combine to consensus (AutoCycler combine)  
9) Polish consensus (minimap2 + racon)  
10) QC assessment (CheckM2)

---

## Requirements

- OS: Linux (recommended). Bash + Conda available in PATH.
- Conda environments:
  - autocycler (activated by the script)
  - checkm2 (activated for QC at the end)

- Tools available inside the autocycler environment:
  - AutoCycler CLI (autocycler)
  - Filtlong (filtlong)
  - Samtools (samtools)
  - Minimap2 (minimap2)
  - Racon (racon)
  - Assemblers you want to use (ensure these are installed and in PATH):
    - canu, flye, metamdbg, miniasm, necat, nextdenovo, plassembler, raven

- Tools available inside the checkm2 environment:
  - CheckM2 (checkm2) + its database file
    - Update the script variable checkm2_database with the absolute path to your CheckM2 database.

Notes:
- CheckM2 requires a separately downloaded database; point checkm2_database to the .dmnd file.

---

## Quick start

1) Put your reads in reads/  
   - Accepted formats: .fastq, .fq, .fastq.gz, .fq.gz, or .bam

2) Edit configurable parameters at the top of the script:
   - THREADS, READ_TYPE, READS_DIR, checkm2_database, FALLBACK_GENOME_SIZE

3) Run:
```bash
chmod +x run_autocycler.sh
./run_autocycler.sh | tee run.log
```

You can re-run the script anytime; it automatically skips completed steps. If a process is interrupted, simply delete the isolate's subdirectory within the incomplete step to force a rerun.

---

## Configuration

Top-of-script parameters you can change:

- THREADS
  - Number of CPU threads to use throughout the pipeline.
- READ_TYPE
  - Set to one of: ont_r9, ont_r10, pacbio_clr, pacbio_hifi
  - Controls minimap2 preset selection for polishing and guides some assembler helpers.
- READS_DIR
  - Folder containing your input read files (default: reads).
- FALLBACK_GENOME_SIZE
  - Used if automatic genome size detection fails (default: 5,000,000).
- checkm2_database
  - Path to your CheckM2 database file (e.g., uniref100.KO.1.dmnd).

Assemblers used
- The script runs (by default) for each sample and each of four subsamples: canu, flye, metamdbg, miniasm, necat, nextdenovo, plassembler, raven
- To adjust, edit the assembler list and/or subsample indices inside Step 2.

---

## Output structure

The script creates a step-wise directory layout with per-isolate subfolders.

- 00_filtered_reads/<sample>/
  - <sample>_converted.fastq (if input was BAM)
  - <sample>_trimmed.fastq (filtlong output)
  - <sample>_genome_size.txt (cached result)
- 01_subsampled_reads/<sample>/
  - sample_01.fastq ... sample_04.fastq
- 02_assemblies/<sample>/
  - {assembler}_{01..04}/* (assembler outputs)
- 03_autocycler_out/<sample>/
  - compressed/ (from autocycler compress)
  - clustering/ (from autocycler cluster; includes qc_pass/cluster_*)
  - consensus_assembly.fasta (from autocycler combine)
- 04_racon_polish/<sample>/
  - genome.racon.fasta (Racon-polished consensus)
- 05_final_consensus/
  - <sample>.fasta (final polished consensus per isolate)
  - checkM2report/quality_report.tsv (CheckM2 outputs)

---

## Resumability and caching

- The script checks for existing files and directories before running each step.
- You can stop and re-run; finished steps are skipped.
- Genome size is cached per isolate in 00_filtered_reads/<sample>_genome_size.txt.

---


## Citations and acknowledgements

Please cite the tools you use in your analyses. This script orchestrates:

- AutoCycler
- Filtlong
- SAMtools
- minimap2
- Racon
- CheckM2
- Assemblers: Canu, Flye, MetaMDBG, miniasm, NECAT, NextDenovo, Plassembler, Raven

Thanks to the authors and maintainers of these excellent tools.

---

## License

This repository contains an orchestration script that wraps external tools. Use of each external tool is subject to its own license. Add your preferred license here for the script itself (e.g., MIT).
