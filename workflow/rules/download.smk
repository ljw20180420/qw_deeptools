rule download_hg38:
    output:
        genome
    log: "results/log/{rules.current.name}.log"
    conda: "../envs/curl.yaml"
    shell:
        "curl --continue-at - --stderr {log} --output {output} {config[genome_url]}"

rule download_rmsk:
    output:
        rmsk
    log: "results/log/{rules.current.name}.log"
    conda: "../envs/curl.yaml"
    shell:
        "curl --continue-at - --stderr {log} --output {output} {config[rmsk_url]}"

rule rmsk_to_bed:
    input:
        rmsk
    output:
        rmsk_bed
    log: "results/log/{rules.current.name}.log"
    conda: "../envs/bedops.yaml"
    shell:
        "zcat {input} | cut -f 2- | tr '\t' ' ' | rmsk2bed | bedops --merge - > {output} 2> {log}"
