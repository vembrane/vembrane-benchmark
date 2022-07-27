rule annotate_variants_snpeff:
    input:
        calls="results/norm/{sample}.vcf",  # .vcf, .vcf.gz or .bcf
        db="resources/snpeff/GRCh38.86",  # path to reference db downloaded with the snpeff download wrapper
    output:
        calls="results/snpeff/vcf_v/{sample}_annotated_snpeff.vcf",  # .vcf, .vcf.gz or .bcf
        stats="results/snpeff/stats/{sample}_variants_snpeff.html",  # summary statistics (in HTML), optional
        csvstats="results/snpeff/stats/{sample}_variants_snpeff.csv",  # summary statistics in CSV, optional
    log:
        "logs/snpeff/{sample}.log",
    params:
        extra="-Xmx4g",  # optional parameters (e.g., max memory 4g)
    wrapper:
        "v1.7.0/bio/snpeff/annotate"


rule snpeff_download:
    output:
        directory("resources/snpeff/GRCh38.86"),
    log:
        "logs/snpeff/download/GRCh38.86.log",
    params:
        reference="GRCh38.86",
    wrapper:
        "v1.7.0/bio/snpeff/download"


rule annotate_variants_vep:
    input:
        calls="results/norm/{sample}.vcf",  # .vcf, .vcf.gz or .bcf
        cache="resources/vep/cache",
        plugins="resources/vep/plugins",
    output:
        calls="results/vep/vcf_v/{sample}_annotated_vep.vcf",  # .vcf, .vcf.gz or .bcf
        stats="results/vep/stats/{sample}_variants_vep.html",
    params:
        extra="--vcf_info_field ANN --everything",
        plugins=["LoFtool"],
    log:
        "logs/vep/annotate_{sample}.log",
    threads: 4
    wrapper:
        "v1.7.0/bio/vep/annotate"


rule download_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    params:
        release=101,
    log:
        "logs/vep/download_vep_plugins.log",
    wrapper:
        "v1.7.0/bio/vep/plugins"


rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species="homo_sapiens",
        build="GRCh38",
        release="101",
    log:
        "logs/vep/cache.log",
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "v1.7.0/bio/vep/cache"
