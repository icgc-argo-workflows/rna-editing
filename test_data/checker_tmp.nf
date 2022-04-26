#!/bin/bash nextflow
/*
 * MIT License
 * 
 * Copyright (c) 2022, ICGC ARGO RNA Editing Group
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * Yang Yang <yang.yang@bjmu.edu.cn>
 * Center for Cancer Bioinformatics, Peking University Cancer Hospital & Institute
 */


/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = 'latest'  // package version

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/rna-editing'
]
default_container_registry = 'ghcr.io'
/********************************************************************/

import java.nio.file.Files
import java.nio.file.Paths
import groovy.json.JsonSlurper

parameter_json = file(params.in)
new groovy.json.JsonSlurper().parseText(parameter_json.text).each { k, v -> params[k] = v }

nextflow.enable.dsl = 2

// universal params go here
params.cpus = 1
params.mem = 32  // GB

params.container_registry = ""
params.container_version = ""
params.container = ""
// dir for outputs, must be set when running in local mode

process samtools2vcf {

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    val sample
    file bam_file
    file ref_genome_fa
    path publish_dir

  output:  // output, make update as needed
    path "${sample}.vcf", emit: vcf

  script:
    """
    samtools mpileup -Q 20 -uf ${ref_genome_fa} ${bam_file} | bcftools call -mv > ${sample}.vcf
    """
}


process samtools_filtering {
  
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    val sample
    file vcf
    val cohort_somatic_maf
    val snp_position_file
    val splicesites_pos_file
    val rmsk_region_file
    val TCGA_somatic_file
    val publish_dir

  output:  // output, make update as needed
    path "${sample}_filtered.vcf", emit: filtered_vcf

  script:
    """
    Rscript ${params.src_dir}/rna_editing/samtools/samtools_filtering.R ${vcf} ${sample} ${publish_dir} ${cohort_somatic_maf} ${snp_position_file} ${splicesites_pos_file} ${rmsk_region_file} ${TCGA_somatic_file}
    """
}

process samtools_filtered_anno {

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    val sample
    file filtered_vcf
    val publish_dir

  output:  // output, make update as needed
    path("${sample}_filtered_anno.vcf"), emit: vcf

  script:
    """
    Rscript ${params.src_dir}/rna_editing/samtools/samtools_anno_summary.R ${filtered_vcf} ${sample} ${params.publish_dir} {REDIportal_file} {params.rmsk_region_file}
    """
}

workflow icgcArgoRnaEditingSamtools{
    samtools2vcf(
      params.sample,
      file(params.bam_file),
      file(params.ref_genome_fa),
      params.publish_dir
    )
    samtools_filtering(
      params.sample,
      samtools2vcf.out.vcf,
      file(params.cohort_somatic_maf),
      file(params.snp_position_file),
      file(params.splicesites_pos_file),
      file(params.rmsk_region_file),
      file(params.TCGA_somatic_file),
      params.publish_dir
    )
    samtools_filtered_anno(
      
    )
}

workflow {
  icgcArgoRnaEditingSamtools()
}

