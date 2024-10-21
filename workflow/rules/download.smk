rule download_genome:
    output:
        get_genome()
    log: "results/log/download_genome.log"
    conda: "../envs/curl.yaml"
    shell:
        "curl --continue-at - --stderr {log} --output {output} {config[genome_url]}"

rule download_rmsk:
    output:
        get_rmsk()
    log: "results/log/download_rmsk.log"
    conda: "../envs/curl.yaml"
    shell:
        "curl --continue-at - --stderr {log} --output {output} {config[rmsk_url]}"

rule rmsk_to_bed:
    input:
        get_rmsk()
    output:
        "results/rmsk.bed"
    log: "results/log/rmsk_to_bed.log"
    conda: "../envs/bedops.yaml"
    shell:
        "zcat {input} | cut -f 2- | tr '\t' ' ' | rmsk2bed | bedops --merge - > {output} 2> {log}"
