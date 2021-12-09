###bam to pileup file; filter
#         {params.samtools} mpileup -Q 20 -uf {input.ref} {input.bam} | awk -F'\\t' '$4>0'| {params.bcftools} call -mv > {output}        
rule samtools2vcf:
    input:
        ref = config['reference']['fasta'][ORGANISM],
        bam = "{prefix}/bam/{sampleid}/{sampleid}_Aligned.sortedByCoord.out.bam",
    threads: 1
    params:
        samtools = config['tools']['samtools'],
        bcftools = config['tools']['bcftools'],
    output:
        '{prefix}/samtools/{sampleid}/{sampleid}.vcf'
    shell:
        """
        {params.samtools} mpileup -Q 20 -uf {input.ref} {input.bam} | {params.bcftools} call -mv > {output}        
        """

rule samtools_filtering:
    input:
        vcf = '{prefix}/samtools/{sampleid}/{sampleid}.vcf',
    threads: 1
    params:
        samtools = config['tools']['samtools'],
        src = config['sourcedir'] + config['tools']['samtools_filtering'],
        outpath = "{prefix}/samtools/{sampleid}",
        somatic_maf = config['resource']['somatic_maf']
    output:
        '{prefix}/samtools/{sampleid}/{sampleid}_filtered.vcf'
    shell:
        """
        Rscript {params.src} {input.vcf} {wildcards.sampleid} {params.outpath} {params.somatic_maf}
        """

rule samtools_filtered_anno:
    input:
        vcf = '{prefix}/samtools/{sampleid}/{sampleid}_filtered.vcf',
    threads: 1
    params:
        samtools = config['tools']['samtools'],
        src = config['sourcedir'] + config['tools']['samtools_filtered_anno'],
        outpath = "{prefix}/samtools/{sampleid}",
        REDIportal_file = "/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/rediportal/known_rnaediting_site.txt",
        rmsk_region_file = "/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/rmsk/rmsk_region_v2.gtf",
    output:
        '{prefix}/samtools/{sampleid}/{sampleid}_filtered_anno.vcf'
    shell:
        """
        Rscript {params.src} {input.vcf} {wildcards.sampleid} {params.outpath} {params.REDIportal_file} {params.rmsk_region_file}
        """

