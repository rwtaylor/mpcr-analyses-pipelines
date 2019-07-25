#!/usr/bin/env nextflow

/* Multiplex pcr calling workflow
  Started September 2017

  @Authors
  Ryan Taylor <ryan@ryantaylor.net>

*/

// Add fastq files to the fastqFiles channel
fastqFiles = Channel.fromFilePairs("${params.fastq_files_path}")
//fastqFiles = Channel.fromPath("${params.fastq_files_path}")
fastqFiles = fastqFiles.view()

// Extract FileIDs from file names. This uses a regex defined in the nextflow.config file
fastqFiles = fastqFiles.map{ sample, reads ->
  fileID = sample =~ params.regex_fastq_ID
  [fileID[0][1], reads[0], reads[1]]
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


fastqFiles = fastqFiles.view()

////////////////////////////////////////////////////////////////////////////////////////////////////
/* Trimming */

// Trim fastq files if the config parameter 'trim_fastq_files' is true
if (params.trim_fastq_files) {

  process Trim_galore {
    publishDir path:"${params.publish_directory}/trimmed-fastqs", mode: "copy", overwrite: true

    tag "${params.output_prefix}-${sampleID}"

    cpus 1
    memory { task.cpus * 4.GB }
    
    input:
    set sampleID, file(r1), file(r2) from fastqFiles

    output:
    set sampleID, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into trimmedFastqs

    """
    /usr/local/bin/trim_galore --paired --length 5 --gzip ${r1} ${r2}
    """
  }

} else {

    trimmedFastqs = fastqFiles

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/* Mapping */

if (params.bam_files_path == null) {
  process Mapping_bwa {
    publishDir path:"${params.publish_directory}/bams", mode: "copy", overwrite: true
    tag "${params.output_prefix}-${sampleID}"

    cpus 2
    memory { task.cpus * 4.GB }

    input:
    set sampleID, file(fq1), file(fq2) from trimmedFastqs

    output:
    file("*.bam") into bwaMappedBams
    file("*.bam.bai") into bamIndexes

    script:
    readGroupString="\"@RG\\tID:${sampleID}\\tSM:${sampleID}\\tLB:${sampleID}\\tPL:illumina\""

    """
    set -eo pipefail
    /usr/local/bin/bwa mem -M -R ${readGroupString} -B 3 -t ${task.cpus} ${params.reference} ${fq1} ${fq2} | \
    /usr/local/bin/samtools view -hu - | /usr/local/bin/samtools sort --threads ${task.cpus} -O bam - > ${sampleID}.bam
    /usr/local/bin/samtools index ${sampleID}.bam
    """
  }
} else {
    bwaMappedBams = Channel.fromPath("${params.bam_files_path}/*.bam")
    bamIndexes = Channel.fromPath("${params.bam_files_path}/*.bam.bai")
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
