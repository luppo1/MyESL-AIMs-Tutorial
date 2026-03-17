# MyESL AIMs Tutorial (YRI vs CEU)

## Overview
This repository provides a step-by-step tutorial for using MyESL (Evolutionary Sparse Learning) to identify Ancestry Informative Markers (AIMs).

We use two highly differentiated populations from the 1000 Genomes Project:
- YRI (Yoruba in Ibadan, Nigeria)
- CEU (Utah residents with European ancestry)

This tutorial allows full reproduction of the analysis using a small example dataset.

---

## Directory Structure

```
MyESL-AIMs-Tutorial/
│
├── data/
│   ├── yri_vs_ceu.vcf.gz
│   ├── yri_vs_ceu.vcf.gz.csi
│   └── 1000G_SampleListWithLocations.txt
│
├── scripts/
│   └── vcf2fasta.py
│
├── MyESL/                # cloned during Step 1
│
├── yri_vs_ceu/           # generated during pipeline
│   ├── fasta/
│   │   ├── *.fasta
│   │   └── alignment_list.txt
│   ├── positions/
│   │   └── *_positions.txt
│
├── results_yri_ceu/
│   └── PSS_classes_summary.txt
│
├── final_AIMs.txt
├── PSS_classes_matched.txt
│
└── README.md
```

---

## Quick Start

```bash
# Clone this repository
git clone https://github.com/luppo1/MyESL-AIMs-Tutorial.git
cd MyESL-AIMs-Tutorial
```

---

## Requirements
- Python 3
- bcftools (must be installed and available in your PATH)

---

## Step 1: Install MyESL

```bash
git clone https://github.com/kumarlabgit/MyESL.git
cd MyESL/
cd bin/
chmod 777 *
cd ..

# Test installation
python3 MyESL.py
```

Expected output:
```
usage: MyESL.py [-h] [--classes CLASSES] [--tree TREE] ...
```

---

## Step 2: Input Data

Example dataset:
```
data/yri_vs_ceu.vcf.gz
```

---

## Step 3: Convert VCF to FASTA

```bash
# Make sure bcftools is already added to PATH
# Example: export PATH=/path/to/bcftools/bin:$PATH

# In my case:
export PATH=/home/tup97263/work/tools/bcftools-1.17/bin:$PATH

# Extract genomic range
start=$(bcftools view ../data/yri_vs_ceu.vcf.gz | grep -v '^#' | head -n 1 | cut -f 2)
end=$(bcftools view ../data/yri_vs_ceu.vcf.gz | grep -v '^#' | tail -n 1 | cut -f 2)

# Convert to FASTA
# Set the chromosome number
chr=22
python3 ../scripts/vcf2fasta.py ../data/yri_vs_ceu.vcf.gz ${chr} --start "$start" --end "$end"

# This creates a directory named yri_vs_ceu containing the fasta and positions files

# Organize outputs
cd yri_vs_ceu
mkdir fasta positions
mv *_positions.txt positions/
mv *.fasta fasta/
```

---

## Step 4: Preprocess FASTA

```bash
cd fasta

# Clean the headers
for file in *.fasta; do
    sed 's/>\([^:]*\):.*/>\1/g' "$file" > tmp && mv tmp "$file"
done

for file in *.fasta; do
    awk '/^>/{print; next} {print substr($0,1,5000)}' "$file" > tmp && mv tmp "$file"
done
```

---

## Step 5: Create Alignment List

```bash
ls "$PWD"/*.fasta > alignment_list.txt

# This contains the paths to the fasta files
```

---

## Step 6: Create Classes File

```bash
# Go back to the main MyESL directory
cd ../../

awk 'BEGIN {OFS="\t"} $2=="CEU"{print $1,-1} $2=="YRI"{print $1,1}' \
../data/1000G_SampleListWithLocations.txt > classes.txt

# This contains a user-defined hypothesis. It has two columns, which are tab-separated. The first column contains species names, and the second column contains the response value for the species (+1/-1).
```

---

## Step 7: Run a basic ESL Model

```bash
python3 MyESL.py yri_vs_ceu/fasta/alignment_list.txt \
    --classes classes.txt \
    --stats_out PGHS \
    --lambda1 0.1 \
    --output results_yri_ceu
```

Parameters:
```
--classes: the user-defined hypothesis file
--stats_out: Output statistics
--lambda1: The site sparsity parameter that ranges from 0 to 1
--output: The name of the output directory where all results from MyESL analysis will be stored

```

---

## Step 8: Map ESL's Output to the Positions obtained from the vcf file

```bash
# ESL's out is located in: MyESL/results_yri_ceu, and the file of interest is PSS_classes_summary.txt 

awk 'FNR==NR{pos[$1]=$2; next} $1 in pos {
    split(pos[$1],a,":");
    print a[1], a[2], $2
}' <(cat yri_vs_ceu/positions/*positions.txt) \
results_yri_ceu/PSS_classes_summary.txt \
> PSS_classes_matched.txt

wc -l PSS_classes_matched.txt results_yri_ceu/PSS_classes_summary.txt

# The number of rows from two files must be the same, i.e. 2990
```

---

## Step 9: Final AIMs

```bash
# Map ESL's output to the vcf file and extract the SNP IDs, Position, Chromosome and add the ESL's PSS scores

awk 'BEGIN{print "Chrom\tSNPID\tPosition\tPSS"}
NR==FNR{key[$1"_"$2]=$3; next}
($1"_"$2) in key {
    print $1, $3, $2, key[$1"_"$2]
}' PSS_classes_matched.txt \
<(bcftools view -H ../data/yri_vs_ceu.vcf.gz) \
> final_AIMs.txt
```

---

## Output

```
final_AIMs.txt
```

---

## Showing the head of final_AIMs.txt

```
Chrom   SNPID        Position    PSS
22      rs226522     45575468    0.000885789020685479
22      rs6007494    45581020    0.0010501827928237617
22      rs114692025  45581870    0.0016331698861904442
22      rs6006951    45583107    0.002064603613689542
22      rs146738389  45583670    0.0010501827928237617
22      rs141850672  45584544    0.002064603613689542
22      rs226517     45584586    0.000885789020685479
22      rs117294954  45584931    0.002064603613689542
22      rs226520     45587746    0.000885789020685479
```

---




## Citation
Doughan et al. (2026) — manuscript in preparation

---

## Contact
albert.doughan@temple.edu
