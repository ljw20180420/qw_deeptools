rule computeGCBias:
    input:
        sorted_bam = "resources/data/{sample}.bam.sorted",
        genome = "resources/hg38.2bit",
        blackListFileName = "results/rmsk.bed"
    params:
        effectiveGenomeSize = 2_913_022_398,
        sampleSize = 50_000_000,
        regionSize = 300
    output:
        GCbiasFrequenciesFile = "results/{sample}.GCfreq.txt",
        biasPlot = "results/{sample}.GCfreq.pdf"
    threads: 2
    log: "results/log/{sample}.computeGCBias.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "computeGCBias --bamfile {input.sorted_bam} --effectiveGenomeSize {params.effectiveGenomeSize} --genome {input.genome} --GCbiasFrequenciesFile {output.GCbiasFrequenciesFile} --fragmentLength {config[fragmentLength]} --sampleSize {params.sampleSize} --blackListFileName {input.blackListFileName} --numberOfProcessors {threads} --biasPlot {output.biasPlot} --regionSize {params.regionSize} 2> {log}"

rule correctGCBias:
    input:
        sorted_bam = "resources/data/{sample}.bam.sorted",
        genome = "resources/hg38.2bit",
        GCbiasFrequenciesFile = "results/{sample}.GCfreq.txt"
    output:
        correctedFile = "results/{sample}.bam.sorted.corrected"