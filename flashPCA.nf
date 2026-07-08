#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.nPCs = params.nPCs ?: 20
params.outdir = params.outdir ?: "results/FlashPCA2_PCA"

// ------------------------------------------------------------
// QC reference
// ------------------------------------------------------------
process qc_norm_ref {
  tag "$chr"

  input:
    tuple val(chr), path(vcf), path(vcf_tbi)

  output:
    tuple val(chr), path("ref_${chr}.qc.vcf.gz"), path("ref_${chr}.qc.vcf.gz.tbi")

  script:
  """
  module load bcftools

  bcftools norm -m -both -f "${params.ref_fasta}" "$vcf" | \\
    bcftools view -q ${params.maf} -Q 0.99 | \\
    bcftools annotate -x INFO,^GT | \\
    bcftools view -S "${params.qc_ref_list}" --force-samples | \\
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ref_${chr}.qc.vcf.gz

  tabix -f ref_${chr}.qc.vcf.gz
  """
}

// ------------------------------------------------------------
// QC study
// ------------------------------------------------------------
process qc_norm_study {
  tag "$chr"

  input:
    tuple val(chr), path(vcf), path(vcf_tbi)

  output:
    tuple val(chr), path("study_${chr}.qc.vcf.gz"), path("study_${chr}.qc.vcf.gz.tbi")

  script:
  """
  module load bcftools

  bcftools norm -m -both -f "${params.ref_fasta}" "$vcf" | \\
    bcftools view -q ${params.maf} -Q 0.99 | \\
    bcftools annotate -x INFO,^GT | \\
    bcftools view -S "${params.qc_study_list}" --force-samples | \\
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o study_${chr}.qc.vcf.gz

  tabix -f study_${chr}.qc.vcf.gz
  """
}

// ------------------------------------------------------------
// Shared variants only, exact CHROM:POS:REF:ALT
// ------------------------------------------------------------
process intersect_shared_variants {
  tag "$chr"

  input:
    tuple val(chr), path(ref), path(ref_tbi), path(study), path(study_tbi)

  output:
    tuple val(chr),
          path("${chr}.shared_ref.vcf.gz"), path("${chr}.shared_ref.vcf.gz.tbi"),
          path("${chr}.shared_study.vcf.gz"), path("${chr}.shared_study.vcf.gz.tbi")

  script:
  """
  module load bcftools
  set -euo pipefail

  bcftools isec -n=2 -c none -w1 -Oz -o ${chr}.shared_ref.vcf.gz "$ref" "$study"
  bcftools isec -n=2 -c none -w2 -Oz -o ${chr}.shared_study.vcf.gz "$ref" "$study"

  tabix -f ${chr}.shared_ref.vcf.gz
  tabix -f ${chr}.shared_study.vcf.gz
  """
}

// ------------------------------------------------------------
// Final QC + LD pruning on reference
// Produces PLINK bfile for FlashPCA2
// ------------------------------------------------------------
process final_qc_and_prune_ref {
  tag "$chr"

  input:
    tuple val(chr), path(vcf), path(vcf_tbi)

  output:
    tuple val(chr),
          path("${chr}.ref.pruned.bed"),
          path("${chr}.ref.pruned.bim"),
          path("${chr}.ref.pruned.fam"),
          path("${chr}.prune.prune.in")

  script:
  """
  module load StdEnv/2020 plink/1.9b_6.21-x86_64 bcftools

  if [ "${params.header_bed}" = "TRUE" ]; then
    tail -n +2 ${params.lowcomplexity_bed} | cut -f1-3 > lowcomplexity.clean.bed
  else
    cut -f1-3 ${params.lowcomplexity_bed} > lowcomplexity.clean.bed
  fi

  bcftools view -T ^lowcomplexity.clean.bed "$vcf" -Oz -o ${chr}.filtered.vcf.gz

  plink --vcf ${chr}.filtered.vcf.gz \\
        --double-id \\
        --real-ref-alleles \\
        --snps-only \\
        --maf ${params.maf} \\
        --geno ${params.geno} \\
        --biallelic-only strict \\
        --make-bed \\
        --out ${chr}.ref.final_qc

  plink --bfile ${chr}.ref.final_qc \\
        --indep-pairwise ${params.prune_window} ${params.prune_step} ${params.prune_r2} \\
        --out ${chr}.prune

  plink --bfile ${chr}.ref.final_qc \\
        --extract ${chr}.prune.prune.in \\
        --make-bed \\
        --out ${chr}.ref.pruned
  """
}

// ------------------------------------------------------------
// Apply reference-pruned SNPs to study
// ------------------------------------------------------------
process prepare_study_pruned {
  tag "$chr"

  input:
    tuple val(chr), path(study_vcf), path(study_tbi), path(prune_list)

  output:
    tuple val(chr),
          path("${chr}.study.pruned.bed"),
          path("${chr}.study.pruned.bim"),
          path("${chr}.study.pruned.fam")

  script:
  """
  module load StdEnv/2020 plink/1.9b_6.21-x86_64

  plink --vcf ${study_vcf} \\
        --double-id \\
        --real-ref-alleles \\
        --extract ${prune_list} \\
        --make-bed \\
        --out ${chr}.study.pruned
  """
}

// ------------------------------------------------------------
// Merge per-chromosome PLINK files
// ------------------------------------------------------------
process merge_plink_ref {
  publishDir "${params.outdir}/plink_reference", mode: "copy"

  input:
    path beds
    path bims
    path fams

  output:
    path "reference_pruned.bed"
    path "reference_pruned.bim"
    path "reference_pruned.fam"

  script:
  """
  module load StdEnv/2020 plink/1.9b_6.21-x86_64

  ls *.bed | sed 's/.bed\$//' | sort -V > merge_list_all.txt
  head -n 1 merge_list_all.txt > first.txt
  tail -n +2 merge_list_all.txt > merge_list.txt

  first=\$(cat first.txt)

  if [ -s merge_list.txt ]; then
    plink --bfile \$first --merge-list merge_list.txt --make-bed --out reference_pruned
  else
    cp \${first}.bed reference_pruned.bed
    cp \${first}.bim reference_pruned.bim
    cp \${first}.fam reference_pruned.fam
  fi
  """
}

process merge_plink_study {
  publishDir "${params.outdir}/plink_study", mode: "copy"

  input:
    path beds
    path bims
    path fams

  output:
    path "study_pruned.bed"
    path "study_pruned.bim"
    path "study_pruned.fam"

  script:
  """
  module load StdEnv/2020 plink/1.9b_6.21-x86_64

  ls *.bed | sed 's/.bed\$//' | sort -V > merge_list_all.txt
  head -n 1 merge_list_all.txt > first.txt
  tail -n +2 merge_list_all.txt > merge_list.txt

  first=\$(cat first.txt)

  if [ -s merge_list.txt ]; then
    plink --bfile \$first --merge-list merge_list.txt --make-bed --out study_pruned
  else
    cp \${first}.bed study_pruned.bed
    cp \${first}.bim study_pruned.bim
    cp \${first}.fam study_pruned.fam
  fi
  """
}

// ------------------------------------------------------------
// FlashPCA2 reference PCA
// ------------------------------------------------------------
process run_flashpca_reference {
  tag "run_flashpca_reference"
  publishDir "${params.outdir}/flashpca_reference", mode: "copy"

  input:
    path bed
    path bim
    path fam

  output:
    path "reference_flashpca.eigenvectors"
    path "reference_flashpca.eigenvalues"
    path "reference_flashpca.loadings"
    path "reference_flashpca.meansd"

  script:
  """
  echo "FlashPCA2 version:"
  ${params.flashpca_bin} --help | head -5
  
  ${params.flashpca_bin} \
    --bfile reference_pruned \
    --ndim ${params.nPCs} \
    --outvec reference_flashpca.eigenvectors \
    --outval reference_flashpca.eigenvalues \
    --outload reference_flashpca.loadings \
    --outmeansd reference_flashpca.meansd \
    -v
  """
}

// ------------------------------------------------------------
// FlashPCA2 study projection
// ------------------------------------------------------------
process project_study_flashpca {
  tag "project_study_flashpca"
  publishDir "${params.outdir}/flashpca_projection", mode: "copy"

  input:
    path study_bed
    path study_bim
    path study_fam
    path loadings
    path meansd

  output:
    path "study_flashpca.projected.txt"

  script:
  """
  ${params.flashpca_bin} \
    --bfile study_pruned \
    --project \
    --inmeansd ${meansd} \
    --outproj study_flashpca.projected.txt \
    --inload ${loadings} \
    -v
  """
}

// ------------------------------------------------------------
// Plot reference-only PCA
// ------------------------------------------------------------
process run_pca_analysis_reference {
  tag "run_pca_analysis_reference"
  publishDir "${params.outdir}/PCA_plots", mode: "copy"

  input:
    path ref_pca_file

  output:
    path "*"

  script:
  """
  module load r
  mkdir -p plot_PCA

  Rscript PCA_flashpca.R \
    ${ref_pca_file} \
    ${params.meta_file} \
    NONE \
    reference_only \
    ${params.k} \
    ${params.n_pcs} \
    ${params.threshold_N}
  """
}


// ------------------------------------------------------------
// Plot reference PCA + projected study samples
// ------------------------------------------------------------
process run_pca_analysis_projection {
  tag "run_pca_analysis_projection"
  publishDir "${params.outdir}/PCA_plots", mode: "copy"

  input:
    path ref_pca_file
    path study_pca_file

  output:
    path "*"

  script:
  """
  module load r
  mkdir -p plot_PCA

  Rscript PCA_flashpca.R \
    ${ref_pca_file} \
    ${params.meta_file} \
    ${params.qc_study_list} \
    projection \
    ${params.k} \
    ${params.n_pcs} \
    ${params.threshold_N} \
    ${study_pca_file}
  """
}

// ------------------------------------------------------------
// WORKFLOW
// ------------------------------------------------------------
workflow {

  def ref_ch = Channel.fromPath(params.input_reference, checkIfExists: true)
    .map { f ->
      def m = (f.name =~ /chr(\\d+)/)
      if (m) {
        def c = m[0][1] as Integer
        if (c >= 1 && c <= 22) tuple(c.toString(), f, file(f.toString() + '.tbi'))
      }
    }
    .filter { it != null }

  def ref_qc = qc_norm_ref(ref_ch)

  def ref_for_pruning
  def study_shared_for_projection

  if (params.run_projection) {

    if (!params.input_study || !params.qc_study_list) {
      error "When --run_projection true, you must provide --input_study and --qc_study_list"
    }

    def study_ch = Channel.fromPath(params.input_study, checkIfExists: true)
      .map { f ->
        def m = (f.name =~ /chr(\\d+)/)
        if (m) {
          def c = m[0][1] as Integer
          if (c >= 1 && c <= 22) tuple(c.toString(), f, file(f.toString() + '.tbi'))
        }
      }
      .filter { it != null }

    def study_qc = qc_norm_study(study_ch)

    def shared = intersect_shared_variants(ref_qc.join(study_qc))

    ref_for_pruning = shared.map { chr, ref_vcf, ref_tbi, study_vcf, study_tbi ->
      tuple(chr, ref_vcf, ref_tbi)
    }

    study_shared_for_projection = shared.map { chr, ref_vcf, ref_tbi, study_vcf, study_tbi ->
      tuple(chr, study_vcf, study_tbi)
    }

  } else {

    ref_for_pruning = ref_qc
  }

  def ref_pruned = final_qc_and_prune_ref(ref_for_pruning)

  def ref_beds = ref_pruned.map { it[1] }.collect()
  def ref_bims = ref_pruned.map { it[2] }.collect()
  def ref_fams = ref_pruned.map { it[3] }.collect()

  def (ref_bed, ref_bim, ref_fam) = merge_plink_ref(ref_beds, ref_bims, ref_fams)

  def (ref_pcs, ref_evals, ref_loadings, ref_meansd) = run_flashpca_reference(ref_bed, ref_bim, ref_fam)

  if (params.run_projection) {

    def prune_lists = ref_pruned.map { chr, bed, bim, fam, prune_list ->
      tuple(chr, prune_list)
    }

    def study_for_prune = study_shared_for_projection.join(prune_lists)
      .map { chr, study_vcf, study_tbi, prune_list ->
        tuple(chr, study_vcf, study_tbi, prune_list)
      }

    def study_pruned = prepare_study_pruned(study_for_prune)

    def study_beds = study_pruned.map { it[1] }.collect()
    def study_bims = study_pruned.map { it[2] }.collect()
    def study_fams = study_pruned.map { it[3] }.collect()

    def (study_bed, study_bim, study_fam) = merge_plink_study(study_beds, study_bims, study_fams)

    def projected = project_study_flashpca(
      study_bed,
      study_bim,
      study_fam,
      ref_loadings,
      ref_meansd
    )

    run_pca_analysis_projection(ref_pcs, projected)

    } else {

        run_pca_analysis_reference(ref_pcs)
    }
}

