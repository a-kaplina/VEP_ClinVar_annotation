# README

This script processes a VCF file (hg38 assembly), adds annotations (VEP, ClinVar, PharmGKB), filters variants by clinical significance, and generates final tables for downstream analysis.

The following files must be present in the same directory as the script:

* `make_long_table_script.py`
* `pharmgkb_map_clean_CAT_all.txt`

## Usage

Specify the path to the folder containing the hg38 VCF file to annotate:

```bash
bash VEP_Clinvar_annot_script.sh ~/absolute/path/VCF
```

---

# Requirements

## Install VEP

Follow the official installation guide:

```bash
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git pull
git checkout release/115
perl INSTALL.pl
```

Official documentation:

https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html

---

## Check hg38 cache files

Verify that the hg38 cache is present:

```bash
ls ~/.vep/homo_sapiens
```

The directory should contain:

```bash
115_GRCh38
```

If it is missing, download the VEP cache:

```bash
cd $HOME/.vep/homo_sapiens
curl -O https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
tar xzf homo_sapiens_vep_115_GRCh38.tar.gz
mv ~/.vep/homo_sapiens/homo_sapiens/115_GRCh38 ~/.vep/homo_sapiens/
```

After installation, the directory should look like:

```bash
ls $HOME/.vep/homo_sapiens
115_GRCh37 115_GRCh38 homo_sapiens_vep_115_GRCh38.tar.gz
```

---

## Download the GRCh38 primary assembly

```bash
cd $HOME/.vep/homo_sapiens/115_GRCh38
wget ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

Convert it to bgzip format:

```bash
zcat ~/.vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | bgzip -c > ~/.vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz
```

---

## Download ClinVar for hg38

```bash
mkdir ~/clinvar38
cd ~/clinvar38
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
bcftools index clinvar.vcf.gz
```

---

# Workflow Overview

Schematic workflow:
<img width="1111" height="908" alt="ENG_VEP_diagram_git drawio" src="https://github.com/user-attachments/assets/8b2e89e3-e102-4988-899b-8dd74002e035" />



---

# Script Workflow Description

## Input validation

* Checks whether the input folder containing source files is provided
* Creates the following folder structure inside it:

  * `input/` (temporary files)
  * `output/`

---

## Steps 1–2: VCF sorting and normalization

* Sorts the input VCF
* Splits multiallelic variants into separate records
* Normalizes variants into a standardized representation

Creates:

```bash
*_norm.vcf.gz
```

---

## Step 3: Create chromosome mapping file

Creates a mapping file used to remove the `chr` prefix from chromosome names in the VCF.

---

## Step 4: Remove `chr` prefix

Converts chromosome names for ClinVar compatibility:

* `chr1 → 1`
* `chrX → X`

Creates:

```bash
*_nochr.vcf.gz
```

---

## Step 5: Remove intermediate `_norm` files

Deletes temporary normalized files.

---

## Step 6: VEP annotation + ClinVar

Runs Ensembl Variant Effect Predictor (VEP) to add:

* Functional consequences (missense, nonsense, etc.)
* SIFT and PolyPhen predictions
* Frequencies from 1000 Genomes and gnomAD
* HGVS nomenclature
* Gene, transcript, and protein domain information

Then adds ClinVar annotations:

* Clinical significance
* Phenotypes
* Review status

Creates:

```bash
*_annotated.vcf.gz
```

---

## Step 7: Remove temporary `_nochr` file

Deletes the temporary file generated after chromosome renaming.

---

## Step 8: Add PharmGKB annotations

Uses the precompiled file:

```bash
pharmgkb_map_clean_CAT_all.txt
```

This file was generated from:

* `ClinicalVariants.tsv` from the ClinPGX website
* Additional variants identified through PharmCAT annotation (marked as `PharmCAT`)

Adds PharmGKB information including evidence levels:

* 1A
* 1B
* 2A
* 2B
* 3
* 4
* PharmCAT-derived variants

Creates:

```bash
*_annotated_Pharm.vcf.gz
```

Also counts variants by PharmGKB evidence level.

---

## Step 9: Remove intermediate annotated file

Deletes:

```bash
*_annotated.vcf.gz
```

keeping only the PharmGKB-annotated version.

---

## Step 10: Genotype filtering

Removes variants where all samples are homozygous reference.

Creates:

```bash
*_real_annotated_Pharm.vcf.gz
```

---

## Step 11: Pathogenicity statistics

Generates:

```bash
*_clnsig_counts.txt
```

These files contain counts of:

* Clinical significance categories
* Review status categories

Statistics are generated both:

* Before filtering
* After removing homozygous reference variants

---

## Step 12: Extract pathogenic variants

Keeps only variants labeled as:

* `Pathogenic`
* `Likely_pathogenic`
* Their combinations

Creates:

```bash
*_all_pathogenic.vcf.gz
```

---

## Step 13: Filter pathogenic variants by review quality

Keeps only pathogenic variants with reliable ClinVar review status:

* 2 stars — multiple submitters, no conflicts
* 3 stars — expert panel
* 4 stars — clinical practice guideline
* 1 star — criteria provided by at least one submitter

Creates:

```bash
*_good_q_pathogenic.vcf.gz
```

---

## Step 14: Extract VUS variants

Filters annotated variants and keeps only:

* `Uncertain_significance`

Creates:

```bash
*_all_vus.vcf.gz
```

---

## Step 15: Filter VUS variants by review quality

Keeps only VUS variants with good review status.

Creates:

```bash
*_good_q_vus.vcf.gz
```

---

## Steps 16–17: Generate tables using `make_long_table_script.py`

Runs the Python script to convert VCF files into tabular format.

Creates:

* `*_Result_pat_fin_annotation.tsv` — pathogenic variants
* `*_Result_vus_fin_annotation.tsv` — VUS variants

The resulting tables contain all annotations, including:

* Population frequencies
* Functional predictions
* Phenotypes
* Drug-related annotations
* PharmGKB information
* Additional metadata

