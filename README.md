# FlashPCA2_AncestryInference

**FlashPCA2_AncestryInference** is a modular [Nextflow DSL2](https://www.nextflow.io/) pipeline for **population structure analysis and ancestry inference** using **FlashPCA2**. The pipeline performs rigorous variant quality control, LD pruning, principal component analysis (PCA) on a reference panel, and optional projection of study samples onto the reference PCA space.

Inspired by Peyton McClelland's pipeline: https://github.com/CERC-Genomic-Medicine/HGDP_1KG_ancestry_inference

---

# 🧬 Pipeline overview

The pipeline supports two execution modes:

### **1. Reference PCA only**

```
Reference VCF
      │
      ▼
 QC & normalization
      │
      ▼
 Remove low-complexity regions
      │
      ▼
 LD pruning
      │
      ▼
 FlashPCA2
      │
      ▼
 PCA plots
```

### **2. Reference PCA + Study Projection**

```
Reference VCF        Study VCF
      │                  │
      ▼                  ▼
 QC & normalization   QC & normalization
      │                  │
      └───────┬──────────┘
              ▼
     Shared variant selection
              │
              ▼
   QC & LD pruning (Reference)
              │
              ▼
       FlashPCA2 Reference PCA
              │
              ▼
 Apply identical SNP set to Study
              │
              ▼
      FlashPCA2 Projection
              │
              ▼
     PCA visualization + KNN ancestry inference
```

---

# Pipeline steps

## 1. Reference VCF normalization and filtering

The reference VCF is normalized using **bcftools**:

- Split multiallelic variants
- Normalize REF/ALT alleles
- Remove INFO fields
- Keep variants passing MAF filters
- Restrict to samples listed in `qc_ref_list`

---

## 2. Study VCF normalization (optional)

When `run_projection = true`, the study VCF undergoes the same normalization procedure before retaining only variants shared with the reference dataset.

---

## 3. Shared variant selection (projection mode only)

Reference and study datasets are intersected using

```
bcftools isec
```

to ensure that only identical variants (same chromosome, position, REF and ALT alleles) are used for projection.

---

## 4. Variant QC and LD pruning

For the reference panel:

- remove low-complexity and MHC regions
- PLINK quality control
- MAF filtering
- missingness filtering
- retain only bi-allelic SNPs
- LD pruning

The resulting SNP set defines the PCA space.

When projection is enabled, the exact same pruned SNP set is extracted from the study dataset.

---

## 5. FlashPCA2

FlashPCA2 computes the PCA from the reference panel and outputs:

- eigenvectors
- eigenvalues
- SNP loadings
- SNP means and standard deviations

When projection is enabled, FlashPCA2 projects the study samples onto the reference PCA using the previously computed SNP loadings.

---

## 6. PCA visualization

The supplied `bin/PCA.R` script:

- plots the reference PCA
- overlays projected study samples (optional)
- performs k-nearest neighbour ancestry assignment
- exports ancestry prediction tables
- produces publication-quality PCA figures

---

# Repository structure

```text
FlashPCA2_AncestryInference/
│
├── flashPCA.nf
├── nextflow.config
├── launch_da.sh
├── README.md
│
├── bin/
│   ├── PCA.R
│   ├── genome_gap_hg38_and_MHC.bed
│
└── results/
```

Scripts inside `bin/` are automatically staged by the pipeline whenever required.

---

# FlashPCA2 installation

FlashPCA2 is not distributed as a standard module on many HPC systems.

Download the latest binary from:

https://github.com/gabraham/flashpca

Specify its location in the configuration file:

```groovy
flashpca_bin = "/path/to/flashpca_x86-64"
```

---

# Inputs

All input files are specified in `nextflow.config`.

| Input | Description |
|--------|-------------|
| `input_reference` | Per-chromosome reference VCFs (`chr*.vcf.gz`) |
| `qc_ref_list` | Reference samples to include |
| `input_study` | Study VCFs (projection mode only) |
| `qc_study_list` | Study sample list (projection mode only) |
| `meta_file` | Sample metadata (ID, population) |
| `ref_fasta` | Uncompressed hg38 FASTA |
| `lowcomplexity_bed` | BED of excluded regions |
| `flashpca_bin` | FlashPCA2 executable |

---

# Outputs

Results are written under

```
results/FlashPCA2_PCA/
```

Typical outputs include

| Output | Description |
|---------|-------------|
| `flashpca_reference/` | FlashPCA2 PCA outputs |
| `flashpca_projection/` | Projected study PCs |
| `plink_reference/` | Pruned reference PLINK files |
| `plink_study/` | Pruned study PLINK files |
| `PCA_plots/` | PCA figures and ancestry assignments |

---

# Parameters

| Parameter | Description | Default |
|------------|-------------|---------|
| `run_projection` | Perform study projection | `false` |
| `nPCs` | Number of PCs | `20` |
| `maf` | Minimum allele frequency | `0.01` |
| `geno` | Maximum missing genotype rate | `0.001` |
| `prune_window` | LD pruning window | `1000` |
| `prune_step` | LD pruning step | `100` |
| `prune_r2` | LD pruning r² threshold | `0.01` |
| `k` | Number of neighbours for KNN ancestry | `10` |
| `n_pcs` | Number of PCs used for KNN | `4` |
| `threshold_N` | Minimum neighbour count for ancestry assignment | `6` |

---

# Dependencies

- Nextflow (DSL2)
- Java ≥17
- bcftools
- PLINK 1.9
- FlashPCA2
- R ≥4.5

Required R packages:

- data.table
- ggplot2
- FNN

---

# Usage

## Reference PCA only

```bash
nextflow run flashPCA.nf \
    -c nextflow.config \
    -resume
```

## Reference PCA + Study projection

Set

```groovy
run_projection = true
```

and provide

```
input_study
qc_study_list
```

Then execute

```bash
nextflow run flashPCA.nf \
    -c nextflow.config \
    -resume
```

---

# Notes

- Input VCFs must be bgzipped and indexed.
- Chromosomes are expected to follow the naming convention `chr1`–`chr22`.
- The reference samples included in the PCA are entirely determined by `qc_ref_list`.
- Intermediate files are cached under `work/` and can be reused with `-resume`.
- Projection uses the identical SNP set selected during reference LD pruning, ensuring reproducible PCA coordinates.

---

# Citation

If you use this pipeline, please cite:

Pelletier J. (2026)

**FlashPCA2_AncestryInference: A Nextflow pipeline for scalable ancestry inference using FlashPCA2.**

GitHub:

https://github.com/JustinPelletier/FlashPCA2_AncestryInference
