# LncAxiom: A Reproducible Snakemake Pipeline for Genome-Wide lncRNA Prediction and Annotation

![Snakefmt](https://img.shields.io/badge/snakefmt-formatted-blue) ![License](https://img.shields.io/badge/license-MIT-green) [![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.8-brightgreen)](https://snakemake.github.io)

**LncAxiom** is a modular, end-to-end Snakemake workflow for the identification, classification, and functional analysis of long non-coding RNAs (lncRNAs) from RNA-seq data. It integrates established tools such as FEELnc, CPC2, HISAT2, StringTie, DESeq2, and psRobot in a reproducible Conda-managed environment.

---
## Workflow DAG

![LncAxiom workflow DAG](DAG.png)

*Overview of the LncAxiom Snakemake workflow as a directed acyclic graph (DAG).*
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
