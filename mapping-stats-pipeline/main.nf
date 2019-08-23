#!/usr/bin/env nextflow

/* Multiplex pcr calling workflow
  Started September 2017

  @Authors
  Ryan Taylor <ryan@ryantaylor.net>

*/

bamFiles = Channel.fromPath("${params.bam_files_path}/*.bam")
baiFiles = Channel.fromPath("${params.bam_files_path}/*.bai")

bamFiles.into{bamFiles; bamFiles1; bamFiles2}
baiFiles.into{baiFiles; baiFiles1; baiFiles2}

process Sample_bam_stats {
  publishDir path:"${params.publish_directory}/sample_bam_stats", mode: "copy", overwrite: true
  tag "${params.output_prefix}"

  cpus 1
  memory { 8.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file(bam) from bamFiles
  file(bais) from baiFiles.toList()

  output:
  file("${bam.baseName}.stats") into bam_stats

  """  
  set -e -o pipefail
  mkdir -p temp
  bamtools stats -in ${bam} > ${bam.baseName}.stats
  """
}
bamFiles1 = bamFiles1.view()

process Taget_coverage {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true
  tag "${params.output_prefix}"

  cpus 1
  memory { 8.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file(bams) from bamFiles1.toList()
  file(bais) from baiFiles1.toList()

  output:
  file("target_sample_coverage.txt") into target_sample_coverage

  script:
  header = bams.collect{"$it.baseName"}.join('\t')
  header = "chrom\tstart\tstop\tsnp\t" + header
  """  
  set -e -o pipefail
  mkdir -p temp
  echo "${header}" > target_sample_coverage.txt
  bedtools multicov -bed ${params.mapping_targets_bed} -bams ${bams} >> target_sample_coverage.txt
  """
}

process Concatenate_bams {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true
  tag "${params.output_prefix}"

  cpus 1
  memory { 8.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file(bams) from bamFiles2.toList()
  file(bais) from baiFiles2.toList()

  output:
  set file("all_samples.bam"), file("all_samples.bam.bai") into all_samples_bam

  script:
  input_bams = bams.collect{"-in $it"}.join(' ')

  """  
  set -e -o pipefail
  mkdir -p temp
  bamtools merge ${input_bams} | bamtools sort -out all_samples.bam
  bamtools index -in all_samples.bam
  """
}

all_samples_bam.into{all_samples_bam; all_samples_bam1; all_samples_bam2; all_samples_bam3}

process Target_coverage_across_samples {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true
  tag "${params.output_prefix}"

  cpus 1
  memory { 8.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set file(allsamplesbam), file(allsamplesbai) from all_samples_bam

  output:
  file("target_coverage.txt") into target_coverage

  """  
  set -e -o pipefail
  mkdir -p temp
  bedtools coverage -a ${params.mapping_targets_bed} -b $allsamplesbam > target_coverage.txt
  """
}

process Ontarget_hits {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true
  tag "${params.output_prefix}"

  cpus 1
  memory { 8.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set file(allsamplesbam), file(allsamplesbai) from all_samples_bam1

  output:
  set file("ontarget_reads.bam"), file("ontarget_reads_stats.txt") into ontarget_reads

  """  
  set -e -o pipefail
  mkdir -p temp
  bedtools intersect -a $allsamplesbam -b ${params.mapping_targets_bed} > ontarget_reads.bam
  bamtools stats -in ontarget_reads.bam > ontarget_reads_stats.txt
  """
}

process Offsite_hits {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true
  tag "${params.output_prefix}"

  cpus 1
  memory { 8.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set file(allsamplesbam), file(allsamplesbai) from all_samples_bam2

  output:
  set file("off_target_reads.bam"), file("off_target_reads_stats.txt") into offtarget_reads

  """  
  set -e -o pipefail
  mkdir -p temp
  bedtools subtract -A -a $allsamplesbam -b ${params.mapping_targets_bed} > off_target_reads.bam
  bamtools stats -in off_target_reads.bam > off_target_reads_stats.txt
  """
}

process Filter_mapping_quality {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true
  tag "${params.output_prefix}"

  cpus 1
  memory { 8.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set file(allsamplesbam), file(allsamplesbai) from all_samples_bam3
  each filter from Channel.from([[">=20","overequal_mq20"],["<20","under_mq20"]])

  output:
  set file("all_samples_${filter[1]}.bam"), file("all_samples_${filter[1]}.targetcoverage.txt") into filtered_bams

  """  
  set -e -o pipefail
  mkdir -p temp
  bamtools filter -mapQuality "${filter[0]}" -in $allsamplesbam -out all_samples_${filter[1]}.bam
  bedtools coverage -a ${params.mapping_targets_bed} -b all_samples_${filter[1]}.bam > all_samples_${filter[1]}.targetcoverage.txt
  """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
