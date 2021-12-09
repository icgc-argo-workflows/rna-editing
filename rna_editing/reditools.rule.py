# detect de-novo editing site
rule reditools_denovo:
    input:
        bam = "{prefix}/bam/{sampleid}/{sampleid}_Aligned.sortedByCoord.out.bam",
        fa = config['reference']['fasta']['hg19']
    output:
        f1 = "{prefix}/reditools_denovo/{sampleid}.denovo.done",
    singularity:
        config['enviroment']['REDItool']
    params:
        outpath = "{prefix}/reditools_denovo",
    shell:
        """
        set +eu
        source activate nature_protocol && \
        python /REDItools/main/REDItoolDenovo.py -t 4 -i {input.bam} -f {input.fa} -o {params.outpath} -F {wildcards.sampleid} && \
        conda deactivate && \
        set -eu && \
        touch {output.f1}
        """

rule reditools_denovo_filtering_inhouse:
    input:
        f1 = "{prefix}/reditools_denovo/{sampleid}.denovo.done",
        rmsk_file = config['resource']['rmsk_file'],
        snp151_file = config['resource']['snp151_file'],
        rediportal_file = config['resource']['rediportal_file'],
        splicesites_gtf_file = config['resource']['splicesites_gtf_file'],
    threads: 1
    params:
        src = config['sourcedir'] + config['tools']['REDItools_dnarna_filtering_inhouse'],
        reditools_path = "{prefix}/reditools_denovo",
        outpath = "{prefix}/reditools_denovo_filtering_inhouse",
        somatic_maf = config['resource']['somatic_maf'],
    output:
        f1 = '{prefix}/reditools_denovo_filtering_inhouse/{sampleid}/{sampleid}_filtering_summary.txt',
        f2 = "{prefix}/reditools_denovo_filtering_inhouse/{sampleid}/outTable_{sampleid}.out.editing.strict",
    singularity:
        config['enviroment']['REDItool']
    shell:
        """
        set +eu
        source activate nature_protocol && \
        sh {params.src} {params.reditools_path} {params.outpath} {wildcards.sampleid} {params.somatic_maf} {input.rmsk_file} {input.snp151_file} {input.rediportal_file} {input.splicesites_gtf_file} && \
        conda deactivate && \
        set -eu
        """

rule reditools_filtering:
    input:
        f1 = "{prefix}/reditools_denovo/{sampleid}.denovo.done",
    threads: 1
    params:
        src = config['sourcedir'] + config['tools']['REDItools_denovo_filtering'],
        reditools_path = "{prefix}/reditools_denovo",
        outpath = "{prefix}/reditools_denovo_filtering"
    output:
        reditools_denovo_file = '{prefix}/reditools_denovo_filtering/{sampleid}/{sampleid}_filtering_summary.txt',
        reditools_denovo_strict_file = "{prefix}/reditools_denovo_filtering/{sampleid}/outTable_{sampleid}.out.editing.strict",
    singularity:
        config['enviroment']['REDItool']
    shell:
        """
        set +eu
        source activate nature_protocol && \
        sh {params.src} {params.reditools_path} {params.outpath} {wildcards.sampleid}
        conda deactivate && \
        set -eu
        """


