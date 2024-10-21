rule index_sorted_bam:
    input:
        "{fullpath}.bam"
    output:
        r"{fullpath}.bam.bai"
    conda: "../envs/bedops.yaml"
    shell:
        "samtools index -b -@ {threads} {input} {output}"

