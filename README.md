# rna-editing

This repository collects the tool packages used for RNA Editing analysis in the context of the ICGC ARGO project. 

# checker.nf
  "container_registry": "ghcr.io",
  "container_version": "latest",
  "container": "ghcr.io/icgc-argo-workflows/rna-editing",
# nextflow run checker_tmp.nf -params-file ./input/SRR388226_chr1.json
# nextflow -C ../nextflow.config run checker_tmp.nf --in ./input/SRR388226_chr1.json
# nextflow -c ../nextflow.config run checker_tmp.nf --in ./input/SRR388226_chr1.json

## succeed
nextflow run checker_tmp.nf --in ./input/SRR388226_chr1.json

nextflow run checker_tmp.nf --in ./input/SRR388226_chr1.json

nextflow run checker.nf --in ./input/SRR388226_chr1.json


# container
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  // publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir


nextflow run checker.nf --in ./input/SRR388226_chr1.json
parameter_json = file(params.in)
new groovy.json.JsonSlurper().parseText(parameter_json.text).each { k, v -> params[k] = v }
include { icgcArgoRnaEditingSamtools } from "../rna_editing/main"


nextflow run checker.nf -params-file ./input/SRR388226_chr1.json


include { icgcArgoRnaSeqAlignmentSTAR } from '../alignSTAR' params(['cleanup': false, *:params])
// include { icgcArgoRnaEditingSamtools } from "../rna_editing/main" params(params)
include { icgcArgoRnaEditingSamtools } from "../rna_editing/main"
// include { icgcArgoRnaEditingRediTools } from "../main" params(params)


workflow icgcArgoRnaEditingSamtools{
  samtools2vcf(
    params.sample,
    file(params.bam_file),
    file(params.ref_genome_fa),
    params.publish_dir
  )
}

  icgcArgoRnaEditingSamtools(
    params.sample,
    file(params.bam_file),
    file(params.ref_genome_fa),
    params.publish_dir
  )
/*
 icgcArgoRnaEditingRediTools(
   params.sequencing_files,
   params.ref_genome_fa,
   params.snp_datasets_file,
   params.cohort_somatic_maf,
   params.rmsk_file,
   params.snp151_file,
   params.rediportal_file,
   params.splicesites_gtf_file
 )
*/


    samtools_filtering(
      file(samtools2vcf.vcf),
      
  )


