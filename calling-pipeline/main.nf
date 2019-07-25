#!/usr/bin/env nextflow

/* Multiplex pcr calling workflow
  Started September 2017

  @Authors
  Ryan Taylor <ryan@ryantaylor.net>

*/

// Add fastq files to the fastqFiles channel
bamFiles = Channel.fromPath("${params.mapped_bams_path}/*.bam")
baiFiles = Channel.fromPath("${params.mapped_bams_path}/*.bai")
//fastqFiles = Channel.fromPath("${params.fastq_files_path}")

/*// Read in regions
regionTasks = Channel
  .from(file(params.regions_file).readLines())
  .map {line ->
    list       = line.split("\t")
    task       = list[0]
    seq_cumsum = list[1]
    flag       = list[2]
    [ task, flag ]
}
*/
targetTasks = Channel
  .from(file(params.snp_locations_tab).readLines())
  .map { line ->
    list     = line.split("\t")
    name     = list[0]
    name     = name.replaceFirst(/:/, "_")
    scaffold = list[1]
    position = list[2]
    region = scaffold + ":" + position
    [ name, region ]
  }

targetTasks = targetTasks.view()

bamFiles.into{ bwaMappedBams; bwaMappedBams_wg}
baiFiles.into{ bamIndexes; bamIndexes_wg}

process Pileup_call_target {
  container = '/zstor/containers/singularity/biobase.img'
  publishDir path:"${params.publish_directory}/vcfs", mode: "copy", overwrite: true
  tag "${params.output_prefix}-${target_name}"

  cpus 1
  memory { 8.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set target_name, target_region from targetTasks
  file(bams) from bwaMappedBams.toList()
  file(bais) from bamIndexes.toList()

  output:
  file("${params.output_prefix}.${target_name}.target.vcf.gz") into target_vcfs
  file("${params.output_prefix}.${target_name}.target.vcf.gz.tbi") into target_vcf_indexes

  """
  set -e -o pipefail
  mkdir -p temp
  /usr/local/bin/bcftools mpileup -r ${target_region} -a INFO/AD,FORMAT/AD,FORMAT/DP -Ou --max-depth 100000 -f ${params.reference} ${bams} |\
   bcftools call -Ou -m | bcftools sort --temp-dir temp -Oz -o ${params.output_prefix}.${target_name}.target.vcf.gz
   tabix -p vcf ${params.output_prefix}.${target_name}.target.vcf.gz
  """
}


target_vcfs = target_vcfs.view()


process ConcatenateTargetVCFs {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true

  cpus {  task.attempt == 1 ? 8: 16  }
  memory { task.attempt == 1 ? 96.GB: 192.GB }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 2
  maxErrors '-1'

  input:
  file(vcfs) from target_vcfs.toList()
  file(vcfindexes) from target_vcf_indexes.toList()

  output:
  set file("${params.output_prefix}.target.vcf.gz" ), file("${params.output_prefix}.target.vcf.gz.tbi") into target_vcf

  script:
  input_vcfs = vcfs.collect{"$it"}.join(' ')

  """
  set -e -o pipefail
  mkdir -p temp
  ls *.target.vcf.gz | awk -F '.' '{print \$1"."\$2}'  | xargs -t -I '{}' cp {}.target.vcf.gz {}.cptarget.vcf.gz
  ls *.target.vcf.gz.tbi | awk -F '.' '{print \$1"."\$2}'  | xargs -t -I '{}' cp {}.target.vcf.gz.tbi {}.cptarget.vcf.gz.tbi
  vcf-concat *.cptarget.vcf.gz | vcf-sort --temporary-directory temp | bgzip -@ ${task.cpus} > ${params.output_prefix}.target.vcf.gz
  tabix -p vcf ${params.output_prefix}.target.vcf.gz
  """
}



/*
process Pileup_call_wg {
  container = '/zstor/containers/singularity/biobase.img'
  publishDir path:"${params.publish_directory}/region_vcfs", mode: "copy", overwrite: true
  tag "${params.output_prefix}"

  cpus 1
  memory { 8.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file(bams) from bwaMappedBams_wg.toList()
  file(bais) from bamIndexes_wg.toList()
  set regionTask, regions from regionTasks

  output:
  set regionTask, file("region_${regionTask}.vcf.bgz"), file("region_*.tbi") into region_wg_VCFs

  """
  set -e -o pipefail
  mkdir -p temp
  /usr/local/bin/bcftools mpileup -Ou ${regions} --max-depth 10000 -f ${params.reference} ${bams} | bcftools call -Ou -m | bcftools sort --temp-dir temp -Oz -o region_${regionTask}.vcf.bgz
  tabix -p vcf region_${regionTask}.vcf.bgz
  """
}

region_wg_VCFs.into { region_wg_VCFs; region_wg_VCFs2 }
region_vcf_files = region_wg_VCFs.map { id, file, fileindex -> file }
region_vcf_indexes = region_wg_VCFs2.map { id, file, fileindex -> fileindex }
all_VCFs = region_vcf_files.toList()
all_indexes = region_vcf_indexes.toList()


process ConcatenateVCFs {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true

  cpus {  task.attempt == 1 ? 8: 16  }
  memory { task.attempt == 1 ? 96.GB: 192.GB }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 2
  maxErrors '-1'

  input:
  file(vcfs) from all_VCFs
  file(vcfindexes) from all_indexes

  output:
  set file("${params.output_prefix}.wg.vcf.gz" ), file("${params.output_prefix}.wg.vcf.gz.tbi") into vcf_wg

  script:
  input_vcfs = vcfs.collect{"$it"}.join(' ')

  """
  set -e -o pipefail
  mkdir -p temp
  ls *.vcf.bgz | awk -F '.' '{print \$1}'  | xargs -t -I '{}' cp {}.vcf.bgz {}.vcf.gz
  ls *.vcf.bgz.tbi | awk -F '.' '{print \$1}'  | xargs -t -I '{}' cp {}.vcf.bgz.tbi {}.vcf.gz.tbi
  vcf-concat region_*.vcf.gz | vcf-sort --temporary-directory temp | bgzip -@ ${task.cpus} > ${params.output_prefix}.wg.vcf.gz
  tabix -p vcf ${params.output_prefix}.wg.vcf.gz
  """
}
vcfs = vcf_wg.mix {vcf_wg; vcfs}
*/

/*process FilterVCF {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true
  container = '/zstor/containers/singularity/biobase.img'
  publishDir "${params.publish_directory}/tsvs", mode: 'copy'

  cpus 1

  input:
  set file(vcf), file(index) from vcf_wg

  output:
  set file("${params.output_prefix}.targets.vcf.gz"), file("${params.output_prefix}.targets.vcf.gz") into filtered_vcf

  """
  cut -f1 ${params.snp_locations_tab} | sed '1d' | sed "s/:/\t/" | bgzip -c > targets.tsv.gz && tabix -s1 -b2 -e2 targets.tsv.gz
  bcftools view -T targets.tsv.gz -O z -o ${params.output_prefix}.targets.vcf.gz $vcf
  tabix -p vcf ${params.output_prefix}.targets.vcf.gz
  """
}
*/


vcfs = target_vcf

process ConvertToTSV {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true
  container = '/zstor/containers/singularity/biobase.img'
  publishDir "${params.publish_directory}/tsvs", mode: 'copy'

  cpus 1

  input:
  set file(vcf), file(index) from vcfs

  output:
  file("*.tsv") into tsv_files

  script:
  output_prefix = vcf.baseName - ~/\.vcf*/

  """
  bcftools query -H -f '%CHROM\t%POS\t%INDEL\t%QUAL\t%REF\t%ALT{0}\t%DP[\t%PL:%GT:%AC]\n' -o ${output_prefix}.tsv $vcf
  """
}




workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
