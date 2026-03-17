# MyESL AIMs Tutorial (YRI vs CEU)

## Overview
This repository provides a step-by-step tutorial for using **MyESL (Evolutionary Sparse Learning)** to identify **Ancestry Informative Markers (AIMs)**.

We demonstrate the pipeline using two highly differentiated populations from the 1000 Genomes Project:
- **YRI** (Yoruba in Ibadan, Nigeria)
- **CEU** (Utah residents with Northern and Western European ancestry)

This tutorial is designed for **reviewers and users** to easily reproduce results using a small example dataset.

---

## Repository Contents

| File | Description |
|------|-------------|
| `vcf2fasta.py` | Script to convert VCF → FASTA |
| `yri_vs_ceu.vcf.gz` | Example merged VCF (subset of SNPs) |
| `yri_vs_ceu.vcf.gz.csi` | Index file for VCF |
| `alignment_list.txt` | List of FASTA files used by MyESL |
| `1000G_SampleListWithLocations.txt` | Sample population labels |
| `classes.txt` | ESL input file with class labels |
| `PSS_classes_summary.txt` | ESL output (positions + scores) |
| `final_AIMs.txt` | Final mapped AIMs (SNP IDs + scores) |

---

## Requirements

- Python 3
- `bcftools`
- MyESL (install below)

---

## Installation

Clone MyESL:

```bash
git clone https://github.com/kumarlabgit/MyESL.git
cd MyESL/