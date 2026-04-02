# LncAxiom: An Automated Computational Pipeline for Genome-Wide lncRNAs Identification, Annotation, and Functional Prioritization

![Snakefmt](https://img.shields.io/badge/snakefmt-formatted-blue) ![License](https://img.shields.io/badge/license-MIT-green) [![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.8-brightgreen)](https://snakemake.github.io)

**LncAxiom** is a modular, end-to-end Snakemake workflow for the identification, classification, and functional analysis of long non-coding RNAs (lncRNAs) from RNA-seq data. It integrates established tools such as FEELnc, CPC2, HISAT2, StringTie, DESeq2, and psRobot in a reproducible Conda-managed environment.

---

## Workflow DAG

![LncAxiom workflow DAG](DAG.png)

*Overview of the LncAxiom Snakemake workflow as a directed acyclic graph (DAG).*

---

## Description

LncAxiom automates the main lncRNA discovery workflow:

1. Data acquisition
2. Preprocessing and QC
3. Alignment
4. Transcript assembly
5. Merging and quantification
6. Candidate filtering
7. Coding potential prediction
8. Classification
9. Differential expression
10. Target analysis
11. Novelty assessment
12. Functional enrichment

The workflow is modular, so selected modules can be run without repeating all upstream steps.

---

## Requirements

- Linux or macOS
- Conda or Mamba
- Snakemake
- Git

---

## Repository structure

- `Snakefile` : main workflow
- `config/` : configuration files
- `rules/` : Snakemake rule files
- `scripts/` : helper scripts
- `envs/` : Conda environment YAML files
- `resources/` : small required resource files
- `test_data/` : example input files

---

## Quick start

```bash
snakemake --use-conda --cores 4
```

---

## Troubleshooting

### 1. FEELnc classifier error: `Can't locate Bio/DB/SeqFeature/Store.pm in @INC`
**Cause:** A required Perl/BioPerl module is missing.  
**Possible fix:** Install the missing Perl module inside the active environment and rerun the affected step.

### 2. CPC2 error: `No executable svm-scale on CPC2 path!`
**Cause:** CPC2 or its `libsvm` dependency is not properly installed or compiled.  
**Possible fix:** Reinstall CPC2 and ensure the required binaries are available in the environment.

### 3. Conda environment creation fails
**Cause:** Missing packages, dependency conflicts, or network issues.  
**Possible fix:** Retry with Mamba, verify channel priority, and check internet connectivity.

### 4. Empty final high-confidence lncRNA output
**Cause:** No transcripts passed the intersection or filtering criteria.  
**Possible fix:** Check the individual FEELnc and CPC2 outputs separately and review filtering thresholds.

### 5. Local FASTQ files are not detected
**Cause:** File names do not exactly match the sample names in `config/samples.tsv`, or `input_method` is not set correctly.  
**Possible fix:** Confirm sample names, file extensions, and `config/config.yaml` settings.

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
