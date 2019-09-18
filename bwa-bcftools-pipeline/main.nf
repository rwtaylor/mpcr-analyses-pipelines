#!/usr/bin/env nextflow

/* Multiplex pcr calling workflow
  Started September 2017

  @Authors
  Ryan Taylor <ryan@ryantaylor.net>

*/

if (params.bam_files_path == null) {

  // Add fastq files to the fastqFiles channel
  fastqFiles = Channel
      .fromFilePairs( params.fastq_files_path, size: -1 )
      .ifEmpty { error "Cannot find any reads matching: ${params.fastq_files_path}" }

  fastqFiles = fastqFiles.view()

  // Extract FileIDs from file names. This uses a regex defined in the nextflow.config file
  fastqFiles = fastqFiles.map{ file_name, reads ->
    fastq_ID = file_name =~ params.regex_fastq_ID
    [fastq_ID[0][1], reads]
  }

  // Update FileID with SampleID if lookup file was given in config
  if(params.sample_name_file != null){
    // Get sample info data
    sample_file = file("${params.sample_name_file}")
    sample_info = Channel.from(sample_file.readLines())
    sample_info = sample_info.map { line ->
      lsplit = line.split('\t')
      [lsplit[0], lsplit[1]]
    }

    // Update FileIDs with SampleIDs
    fastqFiles = fastqFiles.phase(sample_info).map {fq, si ->
      [si[1], fq[1], fq[2]]
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  /* Trimming */

  // Trim fastq files if the config parameter 'trim_fastq_files' is true
  if (params.trim_fastq_files) {

    process Trim_galore {
      publishDir path:"${params.publish_directory}/trimmed-fastqs", mode: "copy", overwrite: true

      tag "${params.output_prefix}-${sampleID}"

      cpus 1
      memory { task.cpus * 4.GB }
    
      input:
      set sampleID, reads from fastqFiles

      output:
      set sampleID, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into trimmedFastqs

      script:

      def single = reads instanceof Path
      if( !single ) {
          """
          trim_galore --paired --length 5 --gzip ${reads}
          """
      }  
      else {
          """
          trim_galore --length 5 --gzip ${reads}
          """
      }

    }

  } else {

      trimmedFastqs = fastqFiles

  }
}

trimmedFastqs = trimmedFastqs.view()

///////////////////////////////////////////////////////////////////////////////
/* Mapping */

if (!params.bam_files_path) {
  print "RUNNING BWA"
  process Mapping_bwa {
    publishDir path:"${params.publish_directory}/bams", mode: "copy", overwrite: true
    tag "${params.output_prefix}-${sampleID}"

    cpus 2
    memory { task.cpus * 4.GB }

    input:
    set sampleID, file(reads) from trimmedFastqs

    output:
    file("*.bam") into mappedBams
    file("*.bam.bai") into bamIndexes

    script:
    readGroupString="\"@RG\\tID:${sampleID}\\tSM:${sampleID}\\tLB:${sampleID}\\tPL:illumina\""

    """
    set -eo pipefail
    bwa mem -M -R ${readGroupString} -B 3 -t ${task.cpus} ${params.reference} ${reads} | samtools sort --threads ${task.cpus} -O bam - > ${sampleID}.bam
    samtools index ${sampleID}.bam
    """
  }
} else {
    mappedBams = Channel.fromPath("${params.bam_files_path}/*.bam")
    bamIndexes = Channel.fromPath("${params.bam_files_path}/*.bam.bai")
}

mappedBams = mappedBams.view()

mappedBams.ifEmpty { error "Bams failure" }

///////////////////////////////////////////////////////////////////////////////
/* Pileup & call */

if (params.targets) {

  target_regions_file = Channel.fromPath(params.targets_file)

  process Pileup_call_target {
    publishDir path:"${params.publish_directory}/vcfs", mode: "copy", overwrite: true
    tag "${params.output_prefix}-${target_name}"

    cpus 2
    memory { 32.GB }
    time { 6.h }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 5
    maxErrors '-1'

    input:
    file(regions_file) from target_regions_file
    file(bams) from mappedBams.toList()
    file(bais) from bamIndexes.toList()

    output:
    set file("${params.output_prefix}.targets.vcf.gz"), file("${params.output_prefix}.targets.vcf.gz.tbi") into target_vcf

    """
    set -e -o pipefail
    mkdir -p temp
    
    bcftools mpileup --regions-file ${regions_file} -a INFO/AD,FORMAT/AD,FORMAT/DP -Ou --max-depth ${params.maximum_depth} -f ${params.reference} ${bams} |\
    bcftools call -Ou -m | bcftools sort --temp-dir temp -Oz -o ${params.output_prefix}.targets.vcf.gz
    
    tabix -p vcf ${params.output_prefix}.targets.vcf.gz
    """
  }

}

if (params.whole_genome && params.regions) {

  regionTasks = Channel
    .from(file(params.regions_file).readLines())
    .map {line ->
      list       = line.split("\t")
      task       = list[0]
      seq_cumsum = list[1]
      flag       = list[2]
      [ task, flag ]
    }
 
  process Pileup_call_regions {
    tag "${params.output_prefix}-${target_name}"

    cpus 2
    memory { 8.GB }
    time { 6.h }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 5
    maxErrors '-1'

    input:
    region from regions
    file(bams) from mappedBams.toList()
    file(bais) from bamIndexes.toList()

    output:
    file("${params.output_prefix}.region.${region[0]}.vcf.gz") into region_vcfs
    file("${params.output_prefix}.region.${region[0]}.vcf.gz.tbi") into region_vcf_indexes

    """
    set -e -o pipefail
    mkdir -p temp
    bcftools mpileup -r ${region[1]} -a INFO/AD,FORMAT/AD,FORMAT/DP -Ou --max-depth ${params.maximum_depth} -f ${params.reference} ${bams} |\
     bcftools +fill-tags call -Ou -m | bcftools sort --temp-dir temp -Oz -o ${params.output_prefix}.region.${region[0]}.vcf.gz
     tabix -p vcf ${params.output_prefix}.region.${region[0]}.vcf.gz
    """
  }

  process ConcatenateVCFs {
    publishDir path:"${params.publish_directory}/vcfs", mode: "copy", overwrite: true

    cpus {  task.attempt == 1 ? 8: 16  }
    memory { task.attempt == 1 ? 96.GB: 192.GB }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 2
    maxErrors '-1'

    input:
    file(vcfs) from region_vcfs
    file(vcfindexes) from region_vcf_indexes

    output:
    set file("${params.output_prefix}.wg.vcf.gz" ), file("${params.output_prefix}.wg.vcf.gz.tbi") into vcf_wg

    script:
    input_vcfs = vcfs.collect{"$it"}.join(' ')

    """
    set -e -o pipefail
    bcftools concat --remove-duplicates -oz --threads ${task.cpus - 1} ${params.output_prefix}.region*.vcf.gz
    tabix -p vcf ${params.output_prefix}.wg.vcf.gz
    """
  }
}

//mkdir temp
//vcf-concat ${params.output_prefix}.region*.vcf.gz | vcf-sort --temporary-directory temp | bgzip -@ ${task.cpus} > ${params.output_prefix}.wg.vcf.gz

  

if (params.whole_genome && !params.regions) {
  process Pileup_call {
    publishDir path:"${params.publish_directory}/vcfs", mode: "copy", overwrite: true
    tag "${params.output_prefix}-${target_name}"

    cpus 2
    memory { 8.GB }
    time { 6.h }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 5
    maxErrors '-1'

    input:
    file(bams) from mappedBams.toList()
    file(bais) from bamIndexes.toList()

    output:
    set file("${params.output_prefix}.wg.vcf.gz" ), file("${params.output_prefix}.wg.vcf.gz.tbi") into vcf_wg

    """
    set -e -o pipefail
    mkdir -p temp
    bcftools mpileup -a INFO/AD,FORMAT/AD,FORMAT/DP -Ou --max-depth ${params.maximum_depth} -f ${params.reference} ${bams} |\
     bcftools +fill-tags call -Ou -m | bcftools sort --temp-dir temp -Oz -o ${params.output_prefix}.wg.vcf.gz
     tabix -p vcf ${params.output_prefix}.wg.vcf.gz
    """
  }
}

if(params.targets && params.whole_genome){
  vcfs = target_vcf.mix(vcf_wg)
}

if(params.targets && !params.whole_genome){
  vcfs = target_vcf
}

if(!params.targets && params.whole_genome){
  vcfs = vcf_wg
}


process ConvertToTSV {
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
  bcftools query -H -f '%CHROM\t%POS\t%INDEL\t%QUAL\t%REF\t%ALT{0}\t%DP[\t%PL:%GT:%DP:%AD]\n' -o ${output_prefix}.tsv $vcf
  """
}



workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}