rule validate_results:
    input:
        md5sums=group_validation_results,
    output:
        f="results/validation/{annotation}/{scenario}/{filetype}/{file}.identical",
    log:
        "logs/validation/check_md5sums/{annotation}.{scenario}.{filetype}.{file}.log",
    conda:
        "../envs/plain.yaml"
    script:
        "../scripts/validate.py"


rule generate_md5sum:
    input:
        "results/filtered/{annotation}/{tool}/{scenario}/{filetype}/{file}.vcf",
    output:
        "results/filtered/{annotation}/{tool}/{scenario}/{filetype}/{file}.sorted.noheader.md5sum",
    conda:
        "../envs/validation.yml"
    log:
        "logs/validation/generate_md5sum/{annotation}.{tool}.{scenario}.{filetype}.{file}.log",
    shell:
        """
        bcftools sort {input} | vembrane table "CHROM, POS, REF, ALT, QUAL" | md5sum > {output} 2> {log}
        """


rule num_written_records:
    input:
        "results/filtered/{annotation}/{tool}/{scenario}/{filetype}/{file}.vcf",
    output:
        "results/counts_out/{annotation}/{tool}/{scenario}/{filetype}/{file}.txt",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/plot/num_written_records/{annotation}.{tool}.{scenario}.{filetype}.{file}.log",
    shell:
        """
        bcftools view -H {input} | wc -l > {output} 2> {log}
        """
