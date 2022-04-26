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

// universal params go here
params.container_registry = ""
params.container_version = ""
params.container = ""
params.src_dir = ""

params.cpus = 1
params.mem = 32  // GB

// tool specific parmas go here, add / change as needed
params.sample = ""
params.bam_file = "NO_FILE1"
params.ref_genome_fa = "NO_FILE2"
// params.snp_datasets_file = "NO_FILE3"
// params.cohort_somatic_maf = "NO_FILE4"
// params.rmsk_file = "NO_FILE5"
// params.snp151_file = "NO_FILE6"
// params.rediportal_file = "NO_FILE7"
// params.splicesites_gtf_file = "NO_FILE8"

// Samtools

process samtools2vcf {

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    val sample
    file bam_file
    file ref_genome_fa
    path publish_dir

  output:  // output, make update as needed
    path '${publish_dir}/${sample}.vcf', emit: vcf

  script:
    """
    samtools mpileup -Q 20 -uf ${ref_genome_fa} ${bam_file} | bcftools call -mv > ${publish_dir}/${sample}.vcf
    """
}


process samtools_filtering {
  
  cpus params.cpus
  memory "${params.mem} GB"

  input:
    val sample
    file vcf
    file cohort_somatic_maf
    path publish_dir

  output:  // output, make update as needed
    path '${publish_dir}/${sample}_filtered.vcf', emit: filtered_vcf

  script:
    """
    Rscript ${params.src_dir}/rna_editing/samtools/samtools_filtering.R ${vcf} ${sample} ${publish_dir} ${cohort_somatic_maf}
    """
}

/*
process samtools_filtered_anno {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    file vcf = '{prefix}/samtools/{sampleid}/{sampleid}.vcf',

  output:  // output, make update as needed
    path("${sample}_filtered_anno.vcf"), emit: vcf

  script:
    """
    Rscript {params.src} {input.vcf} {wildcards.sampleid} {params.outpath} {params.somatic_maf}
    """
}
 
// Reditools
process reditools_denovo {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    file bam = "{prefix}/bam/{sampleid}/{sampleid}_Aligned.sortedByCoord.out.bam",
    file fa = config['reference']['fasta']['hg19']

  output:  // output, make update as needed
    path("{prefix}/reditools_denovo/{sampleid}.denovo.done"), emit: vcf

  script:
    """
    set +eu
    source activate nature_protocol && \
    python /REDItools/main/REDItoolDenovo.py -t 4 -i {input.bam} -f {input.fa} -o {params.outpath} -F {wildcards.sampleid} && \
    conda deactivate && \
    set -eu && \
    touch {output.f1}
    """
}
 
process reditools_denovo_filtering_inhouse {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    file f1 = "{prefix}/reditools_denovo/{sampleid}.denovo.done",
    file rmsk_file = config['resource']['rmsk_file'],
    file snp151_file = config['resource']['snp151_file'],
    file rediportal_file = config['resource']['rediportal_file'],
    file splicesites_gtf_file = config['resource']['splicesites_gtf_file'],

  output:  // output, make update as needed
    path("{prefix}/reditools_denovo_filtering_inhouse/{sampleid}/{sampleid}_filtering_summary.txt"), emit: vcf
    path("{prefix}/reditools_denovo_filtering_inhouse/{sampleid}/outTable_{sampleid}.out.editing.strict"), emit: vcf

  script:
    """
    set +eu
    source activate nature_protocol && \
    sh {params.src} {params.reditools_path} {params.outpath} {wildcards.sampleid} {params.somatic_maf} {input.rmsk_file} {input.snp151_file} {input.rediportal_file} {input.splicesites_gtf_file} && \
    conda deactivate && \
    set -eu
    """
}
 
process result_merge {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:
    file f1 = "{prefix}/reditools_denovo/{sampleid}.denovo.done",
    file rmsk_file = config['resource']['rmsk_file'],
    file snp151_file = config['resource']['snp151_file'],
    file rediportal_file = config['resource']['rediportal_file'],
    file splicesites_gtf_file = config['resource']['splicesites_gtf_file'],

  output:  // output, make update as needed
    path("{prefix}/reditools_denovo_filtering_inhouse/{sampleid}/{sampleid}_filtering_summary.txt"), emit: vcf
    path("{prefix}/reditools_denovo_filtering_inhouse/{sampleid}/outTable_{sampleid}.out.editing.strict"), emit: vcf

  script:
    """
    set +eu
    source activate nature_protocol && \
    sh {params.src} {params.reditools_path} {params.outpath} {wildcards.sampleid} {params.somatic_maf} {input.rmsk_file} {input.snp151_file} {input.rediportal_file} {input.splicesites_gtf_file} && \
    conda deactivate && \
    set -eu
    """
}
*/

// workflow
workflow icgcArgoRnaEditingSamtools{
    samtools2vcf(
      params.sample,
      file(params.bam_file),
      file(params.ref_genome_fa),
      params.publish_dir
    )
    samtools_filtering(
      params.sample,
      file(samtools2vcf.out.vcf),
      file(params.cohort_somatic_maf),
      params.publish_dir
    )
}

/*
workflow icgcArgoRnaEditingReditools{
  samtools2vcf(
    params.sample,
    params.bam_file,
    params.ref_genome_fa,
    params.snp_datasets_file,
    params.cohort_somatic_maf,
    params.rmsk_file,
    params.snp151_file,
    params.rediportal_file,
    params.splicesites_gtf_file
  )
}

workflow icgcArgoRnaEditingMergeResults{
  samtools(
    params.sample,
    params.bam_file,
    params.ref_genome_fa,
    params.snp_datasets_file,
    params.cohort_somatic_maf,
    params.rmsk_file,
    params.snp151_file,
    params.rediportal_file,
    params.splicesites_gtf_file
  )
}
*/



