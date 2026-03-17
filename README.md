# MyESL AIMs Tutorial (YRI vs CEU)

## Overview
This repository provides a step-by-step tutorial for using MyESL (Evolutionary Sparse Learning) to identify Ancestry Informative Markers (AIMs).

We use two highly differentiated populations from the 1000 Genomes Project:
- YRI (Yoruba in Ibadan, Nigeria)
- CEU (Utah residents with European ancestry)

This tutorial allows full reproduction of the analysis using a small example dataset.

---

## Quick Start

```bash
# Clone this repository
git clone https://github.com/luppo1/MyESL-AIMs-Tutorial.git

# Enter the repository
cd MyESL-AIMs-Tutorial
```

---

## Requirements
- Python 3
- bcftools (must be installed and available in your PATH)

---

## Step 1: Install MyESL

```bash
# Clone MyESL repository
git clone https://github.com/kumarlabgit/MyESL.git

# Enter MyESL directory
cd MyESL/

# Move to binary folder
cd bin/

# Give execute permissions
chmod 777 *

# Return
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
cd data

# Extract genomic range
start=$(bcftools view yri_vs_ceu.vcf.gz | grep -v '^#' | head -n 1 | cut -f 2)
end=$(bcftools view yri_vs_ceu.vcf.gz | grep -v '^#' | tail -n 1 | cut -f 2)

# Convert to FASTA
python3 ../scripts/vcf2fasta.py yri_vs_ceu.vcf.gz 22 --start "$start" --end "$end"

# Organize outputs
mkdir fasta positions
mv *_positions.txt positions/
mv *.fasta fasta/
```

---

## Step 4: Preprocess FASTA

```bash
cd fasta

# Clean headers
for file in *.fasta; do
    sed 's/>\([^:]*\):.*/>\1/g' "$file" > tmp && mv tmp "$file"
done
for file in *.fasta; do
    awk '/^>/{print; next} {print substr($0,1,5000)}' "$file" > tmp && mv tmp "$file"
done
```

---

## Step 5: Alignment List

```bash
ls "$PWD"/*.fasta > alignment_list.txt
```

---

## Step 6: Classes File

```bash
awk 'BEGIN {OFS="\t"} $2=="CEU"{print $1,-1} $2=="YRI"{print $1,1}' \
../1000G_SampleListWithLocations.txt > classes.txt
```

---

## Step 7: Run MyESL

```bash
cd ../..

python3 MyESL.py data/fasta/alignment_list.txt \
    --classes data/fasta/classes.txt \
    --stats_out PGHS \
    --lambda1 0.1 \
    --output results_yri_ceu
```

---

## Step 8: Map Positions

```bash
awk 'FNR==NR{pos[$1]=$2; next} $1 in pos {
    split(pos[$1],a,":");
    print a[1], a[2], $2
}' <(cat data/positions/*positions.txt) \
results_yri_ceu/PSS_classes_summary.txt \
> PSS_classes_matched.txt

wc -l PSS_classes_matched.txt results_yri_ceu/PSS_classes_summary.txt
```

---

## Step 9: Final AIMs

```bash
awk 'BEGIN{print "Chrom\tSNPID\tPosition\tPSS"}
NR==FNR{key[$1"_"$2]=$3; next}
($1"_"$2) in key {
    print $1, $3, $2, key[$1"_"$2]
}' PSS_classes_matched.txt \
<(bcftools view -H data/yri_vs_ceu.vcf.gz) \
> final_AIMs.txt
```

---

## Output

```
final_AIMs.txt
```

---

## Validation

```bash
diff final_AIMs.txt results/final_AIMs.txt
```

---

## Citation
Doughan et al.

---

## Contact
Open a GitHub issue for questions.
