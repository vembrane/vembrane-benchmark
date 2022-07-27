
rule convert_to_bcfb:
    input:
        "results/{annotation}/vcf_v/{file}.vcf",
    output:
        compressed_bcf="results/{annotation}/bcf_b/{file}.bcf",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/download/convert/to_bcfb/{annotation}_{file}.log",
    shell:
        """
        bcftools view {input} --compression-level 9 --output-type b > {output.compressed_bcf} 2> {log}
        """


rule convert_to_bcfu:
    input:
        "results/{annotation}/vcf_v/{file}.vcf",
    output:
        uncompressed_bcf="results/{annotation}/bcf_u/{file}.bcf",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/download/convert/to_bcfu/{annotation}_{file}.log",
    shell:
        """
        bcftools view {input} --output-type u > {output.uncompressed_bcf} 2> {log}
        """


rule convert_to_vcfz:
    input:
        "results/{annotation}/vcf_v/{file}.vcf",
    output:
        compressed_vcf="results/{annotation}/vcf_z/{file}.vcf.gz",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/download/convert/to_vcfz/{annotation}_{file}.log",
    shell:
        """
        bcftools view {input} --output-type z > {output.compressed_vcf} 2> {log}
        """


rule count_num_records:
    input:
        calls="results/norm/{sample}.vcf",
    output:
        counts="results/counts/{sample}.txt",
    benchmark:
        "benchmarks/bcftools_view/{sample}.txt"
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/download/count_num_records.{sample}.log",
    shell:
        """
        bcftools view --no-header {input.calls} | wc -l > {output.counts} 2> {log}
        """


rule bcftools_merge:
    input:
        a="results/norm/{a}.vcf.gz",
        a_i="results/norm/{a}.vcf.gz.csi",
        b="results/norm/{b}.vcf.gz",
        b_i="results/norm/{b}.vcf.gz.csi",
    output:
        ab="results/norm/{a}_{b}.vcf",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/download/merge.{a}_{b}.log",
    shell:
        """
        bcftools merge {input.a} {input.b} -m both | bcftools norm -m-any -Ov > {output.ab} 2> {log}
        """


rule index_bgzip:
    input:
        "results/norm/{file}.vcf.gz",
    output:
        "results/norm/{file}.vcf.gz.csi",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/download/bgzip-index_{file}.log",
    shell:
        """
        bcftools index {input} 2> {log}
        """


rule bgzip:
    input:
        "results/norm/{file}.vcf",
    output:
        "results/norm/{file}.vcf.gz",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/download/bgzip_{file}.log",
    shell:
        """
        bgzip {input} -c > {output} 2> {log}
        """


rule bcftools_norm:
    input:
        calls="results/raw/chr1/{file}.vcf",
    output:
        norm="results/norm/{file}.vcf",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/download/bcftools_norm/{file}.log",
    shell:
        """
        bcftools norm -N -m-any {input.calls} > {output.norm} 2> {log}
        """


rule restrict_to_chromosome_1:
    input:
        "results/raw/{file}.vcf",
    output:
        "results/raw/chr1/{file}.vcf",
    conda:
        "../envs/vembrane.yml"
    log:
        "logs/restrict_to_chromosome_1/{file}.log",
    shell:
        """
        vembrane filter "CHROM == 'chr1'" {input} > {output}
        """


rule gunzip_vcfs:
    input:
        "results/vcfgz/{file}.vcf.gz",
    output:
        "results/raw/{file}.vcf",
    conda:
        "../envs/SnpSift.yml"
    log:
        "logs/download/gunzip_{file}.log",
    shell:
        """gzip -dc {input} | SnpSift filter "(!(ALT has '*'))" > {output} 2> {log}"""


rule download_vcfs:
    output:
        hg001="results/vcfgz/HG001.vcf.gz",
        hg002="results/vcfgz/HG002.vcf.gz",
        hg003="results/vcfgz/HG003.vcf.gz",
        hg004="results/vcfgz/HG004.vcf.gz",
    conda:
        "../envs/curl.yaml"
    log:
        "logs/download/giab_vcfs.log",
    shell:
        """
        curl -L ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz > {output.hg001} 2> {log}
        curl -L ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.1_draft_benchmark.vcf.gz > {output.hg002} 2>> {log}
        curl -L ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.vcf.gz > {output.hg003} 2>>{log}
        curl -L ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.vcf.gz > {output.hg004} 2>>{log}
        """
