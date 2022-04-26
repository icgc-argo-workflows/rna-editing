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
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    val sample
    file bam_file
    file ref_genome_fa
    path publish_dir

  output:  // output, make update as needed
    path "${publish_dir}/${sample}.vcf", emit: vcf

  script:
    """
    samtools mpileup -Q 20 -uf ${ref_genome_fa} ${bam_file} | bcftools call -mv > ${publish_dir}/${sample}.vcf 2>&1     
    """
}

workflow icgcArgoRnaEditingSamtools{
  samtools2vcf(
    params.sample,
    file(params.bam_file),
    file(params.ref_genome_fa),
    params.publish_dir
  )
}

workflow {
  icgcArgoRnaEditingSamtools()
}

