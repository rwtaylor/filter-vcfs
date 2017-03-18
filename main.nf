#!/usr/bin/env nextflow

input_vcfs = Channel.from(params.input_vcfs).map{it -> file(it) }.map{ file -> [file.baseName, file] }

process FilterVCF {
  publishDir "${params.publish_dir}/filtered", mode: 'copy'
  tag {"q" + filterset.minqual + "-gq" + filterset.mingq}
  cpus 1
  memory 4.GB
  time 6.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(vcf) from input_vcfs
  each filterset from params.filtersets
  
  output:
  set val("${prefix}-snp-${filterset.prefix}"), file("*.vcf") into filtered_vcfs

  script:
  removeindividuals = params.excludesamples.collect{"--remove-indv $it"}.join(" ")

  """
  /usr/local/bin/vcftools --vcf ${vcf} --recode --recode-INFO-all --remove-indels ${removeindividuals}\
  --minQ ${filterset.minqual} --minGQ ${filterset.mingq} --hwe ${filterset.hwe}\
  --out ${prefix}-snp-${filterset.prefix}
  """
}

process SortVCF {
  tag {vcf}
  cpus 1
  memory 4.GB
  time 6.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(vcf) from filtered_vcfs

  output:
  set val("${prefix}"), file("*.vcf") into sorted_vcfs

  """
  set -e
  mkdir -p temp
  /usr/local/bin/vcf-sort --temporary-directory temp < ${vcf} > ${prefix}-sorted.vcf
  """
}

sorted_vcfs.into{ sorted_vcfs; sorted_vcfs_to_sample }


process SampleVCF {
  publishDir "${params.publish_dir}/subsampled", mode: 'copy'
  tag {subsamplerate}
  cpus 1
  memory 4.GB
  time 6.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(vcf) from sorted_vcfs_to_sample
  each subsamplerate from params.subsamplerates

  output:
  set val("${prefix}-ss${subsamplerate}"), file("*.vcf") into subsampled_vcfs

  """
  /usr/local/opt/vcflib/bin/vcfrandomsample -r ${subsamplerate} ${vcf} > ${prefix}-ss${subsamplerate}.vcf
  """
}

vcfs_to_rename = sorted_vcfs.mix(subsampled_vcfs)

process RenameChromosomes {
  publishDir "${params.publish_dir}/chr-renamed", mode: 'copy'
  cpus 1
  memory 4.GB
  time 6.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(vcf) from vcfs_to_rename
  
  output:
  set val("${prefix}-chrename"), file("*.vcf") into renamed_vcfs

  """
  perl -p -e 's/N[W,C]_([0-9]*)\\./\$1/g' ${vcf} > ${prefix}-chrename.vcf
  """
}

workflow.onComplete {
  println "Pipeline completed at: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
