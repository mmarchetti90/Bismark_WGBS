#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Pipeline to align Whole Genome Bisulfite Sequencing paired-end experiments using Bismark.
This pipeline will output aligned bam files as well as txt files for CpG/CHG/CHH methylation.
N.B. The pipeline assumes the user is starting from raw fastq paired-end files and that the bismark_genome_preparation has already been run (see config file)
*/

// ---------------------------------------- //

process TrimFastQ {

  label 'ucgd'

  //publishDir "${projectDir}/${params.trimmed_fastq_dir}", mode: "copy", pattern: "${read_id}_val_{1,2}.fq.gz"
  publishDir "${projectDir}/${params.reports_dir}", mode: "copy", pattern: "*_fastqc.{html,zip}"
  publishDir "${projectDir}/${params.reports_dir}", mode: "copy", pattern: "*_trimming_report.txt"
  
  input:
  tuple val(read_id), path(read1), path(read2)

  output:
  path "*_fastqc.{html,zip}"
  path "*_trimming_report.txt"
  tuple val("${read_id}"), path("${read_id}_val_1.fq.gz"), path("${read_id}_val_2.fq.gz"), emit: trimmed_fastq_files

  """
  trim_galore \
  --cores 4 \
  --output_dir . \
  --basename ${read_id} \
  --paired \
  --length 20 \
  --fastqc \
  --gzip \
  ${read1} ${read2}
  """

}

// ---------------------------------------- //

process RunBismarkAlignment {

  label 'ucgd_long'

  //publishDir "${projectDir}/${params.bam_outdir}", mode: "copy", pattern: "*.bam"
  publishDir "${projectDir}/${params.reports_dir}", mode: "copy", pattern: "*.txt"

  input:
  tuple val(read_id), path(read1), path(read2)

  output:
  path "*_bismark_bt2_pe.bam", emit: aligned_bam_files
  path "*_bismark_bt2_PE_report.txt", emit: bismark_reports

  """
  ${params.bismark_path}/bismark \
  --parallel 8 \
  --output_dir . \
  --bowtie2 \
  --bam \
  ${params.genome_files} \
  -1 ${read1} \
  -2 ${read2}
  """

}

// ---------------------------------------- //

process BismarkDeduplicate {

  label 'ucgd'

  publishDir "${projectDir}/${params.deduplicated_outdir}", mode: "copy", pattern: "*.bam"
  publishDir "${projectDir}/${params.reports_dir}", mode: "copy", pattern: "*deduplication_report.txt"

  input:
  path bam_file

  output:
  path "*_bismark_bt2*.deduplicated.bam", emit: deduplicated_bam_files
  path "*deduplication_report.txt", emit: deduplication_reports

  """
  ${params.bismark_path}/deduplicate_bismark \
  --output_dir . \
  --bam \
  --paired \
  ${bam_file}
  """

}

// ---------------------------------------- //

process ExtractMethylation {

  label 'ucgd'

  publishDir "${projectDir}/${params.methylation_outdir}", mode: "copy", pattern: "{CpG,CHG,CHH}*.txt.gz"
  publishDir "${projectDir}/${params.methylation_outdir}", mode: "copy", pattern: "*{cov,bedGraph}*"
  publishDir "${projectDir}/${params.reports_dir}", mode: "copy", pattern: "*{deduplicated.M-bias.txt,deduplicated_splitting_report.txt}"

  input:
  path bam_file

  output:
  path "CpG*.txt.gz"
  path "CHG*.txt.gz"
  path "CHH*.txt.gz"
  path "*cov*"
  path "*bedGraph*"
  path "*deduplicated.M-bias.txt", emit: mbias_reports
  path "*deduplicated_splitting_report.txt", emit: splitting_reports

  """
  ${params.bismark_path}/bismark_methylation_extractor \
  --parallel 8 \
  --output_dir . \
  --paired-end \
  --no_overlap \
  --cytosine_report \
  --genome_folder ${params.genome_files} \
  --bedGraph \
  --gzip \
  ${bam_file}
  """

}

// ---------------------------------------- //

process MakeReports {

  label 'ucgd'

  publishDir "${projectDir}/${params.reports_dir}", mode: "copy", pattern: "*.html"

  input:
  path alignment_reports
  path deduplication_reports
  path mbias_reports
  path splitting_reports

  output:
  path "*.html"

  """
  ls *_bismark_bt2_PE_report.txt | xargs -n 1 -I {} ${params.bismark_path}/bismark2report --alignment_report "{}"
  """

}

// ----------------Workflow---------------- //

workflow {

  // CREATING FASTQ CHANNEL --------------- //

  fastq_files = Channel.fromFilePairs(params.fastq_files, flat: true)

  // TRIM GALORE! ------------------------- //

  TrimFastQ(fastq_files)

  // BISMARK ALIGNMENT -------------------- //

  RunBismarkAlignment(TrimFastQ.out.trimmed_fastq_files)

  // DEDUPLICATE -------------------------- //

  BismarkDeduplicate(RunBismarkAlignment.out.aligned_bam_files)

  // EXTRACT METHYLATION ------------------ //

  ExtractMethylation(BismarkDeduplicate.out.deduplicated_bam_files)

  // CREATING HTML REPORTS ---------------- //

  MakeReports(RunBismarkAlignment.out.bismark_reports.collect(), BismarkDeduplicate.out.deduplication_reports.collect(), ExtractMethylation.out.mbias_reports.collect(), ExtractMethylation.out.splitting_reports.collect())

}