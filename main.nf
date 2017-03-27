#!/usr/bin/env nextflow

input_vcfs = Channel.from(params.input_vcfs).map{it -> file(it) }.map{ file -> [file.baseName, file] }

process FilterVCF {
  publishDir "${params.publish_dir}/filtered", mode: 'copy'
  tag { prefix + "-snp-" + filterset.prefix }
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
  tag { prefix }
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
  tag { prefix + "-ss" + subsamplerate }
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

// Plink requires scaffolds to have a character prefix, so if they are integers then rename
// prepending "0" makes the scaff a non-integer (at least as far as plink is concerned), while
// also allowing LSAK to interpret the scaff ID as an integer... tricky but convenient.

process RenameChromosomes {
  tag { prefix }
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
  # 1) Get scaffold IDs, and output scaffolds plus renamed scaffold to text file
  #    This strips all non-numeric characters from the scaffold ID. May cause issues if non-numeric characters
  #    are important for name uniquness...
  set -e
  grep -oP '^##contig=<ID=.*' $vcf | \
    sed -e 's/^##contig=<ID=\\(.*\\),.*/\\1/gm' | \
    awk '{ printf \$1 " "; gsub(/[A-Z_.]/,"", \$1); print 0\$1}' \
    > scaffs.txt
  # 2) Use BCFtools to rename chromosomes in VCF
  /usr/local/bin/bcftools annotate --rename-chrs scaffs.txt -o $prefix-chr-renamed.vcf $vcf
  """
}

workflow.onComplete {
  println "Pipeline completed at: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
