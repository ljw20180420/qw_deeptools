rule plotFingerprint:
    params:
        ignoreDuplicates = config["plotFingerprint"]["ignoreDuplicates"],
        minMappingQuality = config["plotFingerprint"]["minMappingQuality"],
        samFlagExclude = config["plotFingerprint"]["samFlagExclude"],
        samples = get_all_samples,
        binSize = config["plotFingerprint"]["binSize"],
        numberOfSamples = lambda wildcards: get_sampleSize(wildcards, limit=config["plotFingerprint"]["numberOfSamples"]),
        skipZeros = config["plotFingerprint"]["skipZeros"],
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else ""
    input:
        bamfiles = get_all_sorted_bam,
        bamfiles_index = lambda wildcards: get_all_sorted_bam(wildcards, suffix=".sorted.bam.bai"),
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    output:
        plotFile = "results/plotFingerprint.pdf",
        outRawCounts = "results/plotFingerprint.tsv",
        outQualityMetrics = "results/plotFingerprint.outQualityMetrics"
    log: "results/log/plotFingerprint.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "plotFingerprint --bamfiles {input.bamfiles} --plotFile {output.plotFile} --outRawCounts {output.outRawCounts} {params.ignoreDuplicates} --minMappingQuality {params.minMappingQuality} --samFlagExclude {params.samFlagExclude} --labels {params.samples} --binSize {params.binSize} --numberOfSamples {params.numberOfSamples} {params.skipZeros} --outQualityMetrics {output.outQualityMetrics} {params.blackListFileName} {input.blackListFileName} --numberOfProcessors {threads} 2> {log}"

rule bamPEFragmentSize:
    params:
        samples = get_all_samples,
        binSize = config["bamPEFragmentSize"]["binSize"],
        distanceBetweenBins = config["bamPEFragmentSize"]["distanceBetweenBins"],
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else ""
    input:
        bamfiles = get_all_sorted_bam,
        bamfiles_index = lambda wildcards: get_all_sorted_bam(wildcards, suffix=".sorted.bam.bai"),
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    output:
        histogram = "results/bamPEFragmentSize.pdf",
        table = "results/bamPEFragmentSize.tsv",
        outRawFragmentLengths = "results/outRawFragmentLengths.tsv"
    log: "results/log/bamPEFragmentSize.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "bamPEFragmentSize --bamfiles {input.bamfiles} --histogram {output.histogram} --numberOfProcessors {threads} --samplesLabel {params.samples} --binSize {params.binSize} --distanceBetweenBins {params.distanceBetweenBins} {params.blackListFileName} {input.blackListFileName} --table {output.table} --outRawFragmentLengths {output.outRawFragmentLengths} 2> {log}"

rule plotCoverage:
    params:
        samples = get_all_samples,
        skipZeros = config["plotCoverage"]["skipZeros"],
        numberOfSamples = lambda wildcards: get_sampleSize(wildcards, limit=config["plotCoverage"]["numberOfSamples"]),
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else "",
        ignoreDuplicates = config["plotCoverage"]["ignoreDuplicates"],
        minMappingQuality = config["plotCoverage"]["minMappingQuality"],
        samFlagExclude = config["plotCoverage"]["samFlagExclude"]
    input:
        bamfiles = get_all_sorted_bam,
        bamfiles_index = lambda wildcards: get_all_sorted_bam(wildcards, suffix=".sorted.bam.bai"),
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    output:
        plotFile = "results/plotCoverage.pdf",
        outRawCounts = "results/plotCoverage.tsv"
    log: "results/log/plotCoverage.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "plotCoverage --bamfiles {input.bamfiles} --plotFile {output.plotFile} --labels {params.samples} {params.skipZeros} --numberOfSamples {params.numberOfSamples} --outRawCounts {output.outRawCounts} {params.blackListFileName} {input.blackListFileName} --numberOfProcessors {threads} {params.ignoreDuplicates} --minMappingQuality {params.minMappingQuality} --samFlagExclude {params.samFlagExclude} 2> {log}"

rule plotEnrichment:
    params:
        samples = get_all_samples,
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else "",
        ignoreDuplicates = config["plotEnrichment"]["ignoreDuplicates"],
        minMappingQuality = config["plotEnrichment"]["minMappingQuality"],
        samFlagExclude = config["plotEnrichment"]["samFlagExclude"]
    input:
        bamfiles = get_all_sorted_bam,
        bamfiles_index = lambda wildcards: get_all_sorted_bam(wildcards, suffix=".sorted.bam.bai"),
        BED = get_all_peaks,
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    output:
        plotFile = "results/plotEnrichment.pdf",
        outRawCounts = "results/plotEnrichment.tsv"
    log: "results/log/plotEnrichment.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "plotEnrichment --bamfiles {input.bamfiles} --BED {input.BED} --plotFile {output.plotFile} --labels {params.samples} --smartLabels --outRawCounts {output.outRawCounts} {params.blackListFileName} {input.blackListFileName} --numberOfProcessors {threads} {params.ignoreDuplicates} --minMappingQuality {params.minMappingQuality} --samFlagExclude {params.samFlagExclude}"

rule estimateReadFiltering:
    params:
        samples = get_all_samples,
        binSize = config["estimateReadFiltering"]["binSize"],
        distanceBetweenBins = config["estimateReadFiltering"]["distanceBetweenBins"],
        ignoreDuplicates = config["alignmentSieve"]["ignoreDuplicates"], # use the same setting for alignmentSieve
        minMappingQuality = config["alignmentSieve"]["minMappingQuality"],
        samFlagExclude = config["alignmentSieve"]["samFlagExclude"],
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else ""
    input:
        bamfiles = get_all_sorted_bam,
        bamfiles_index = lambda wildcards: get_all_sorted_bam(wildcards, suffix=".sorted.bam.bai"),
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    output:
        outFile = "results/estimateReadFiltering.tsv"
    log: "results/log/estimateReadFiltering.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "estimateReadFiltering --bamfiles {input.bamfiles} --outFile {output.outFile} --sampleLabels {params.samples} --binSize {params.binSize} --distanceBetweenBins {params.distanceBetweenBins} --numberOfProcessors {threads} {params.ignoreDuplicates} --minMappingQuality {params.minMappingQuality} --samFlagExclude {params.samFlagExclude} {params.blackListFileName} {input.blackListFileName} 2> {log}"

rule alignmentSieve:
    params:
        genomeChunkLength = config["alignmentSieve"]["genomeChunkLength"],
        ignoreDuplicates = config["alignmentSieve"]["ignoreDuplicates"],
        minMappingQuality = config["alignmentSieve"]["minMappingQuality"],
        samFlagExclude = config["alignmentSieve"]["samFlagExclude"],
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else ""
    input:
        bam = get_sorted_bam,
        bam_index = lambda wildcards: get_sorted_bam(wildcards, suffix=".sorted.bam.bai"),
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    log: "results/log/alignmentSieve.{sample}.log"
    conda: "../envs/deeptools.yaml"
    output:
        outFile = "results/{sample}.filtered.sorted.bam",
        filterMetrics = "results/{sample}.filterMetrics",
        filteredOutReads = "results/{sample}.filteredOutReads"
    shell:
        "alignmentSieve --bam {input.bam} --outFile {output.outFile} --numberOfProcessors {threads} --filterMetrics {output.filterMetrics} --filteredOutReads {output.filteredOutReads} --label {wildcards.sample} --genomeChunkLength {params.genomeChunkLength} {params.ignoreDuplicates} --minMappingQuality {params.minMappingQuality} --samFlagExclude {params.samFlagExclude} {params.blackListFileName} {input.blackListFileName} 2> {log}"

rule computeGCBias:
    params:
        effectiveGenomeSize = get_effectiveGenomeSize,
        fragmentLength = config["fragmentLength"],
        sampleSize = lambda wildcards: get_sampleSize(wildcards, limit=config["computeGCBias"]["sampleSize"]),
        regionSize = config["computeGCBias"]["regionSize"],
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else ""
    input:
        bamfile = "results/{sample}.filtered.sorted.bam",
        bamfile_index = "results/{sample}.filtered.sorted.bam.bai",
        genome = get_genome(),
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    output:
        GCbiasFrequenciesFile = "results/{sample}.GCfreq.txt",
        biasPlot = "results/{sample}.GCfreq.pdf"
    log: "results/log/computeGCBias.{sample}.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "computeGCBias --bamfile {input.bamfile} --effectiveGenomeSize {params.effectiveGenomeSize} --genome {input.genome} --GCbiasFrequenciesFile {output.GCbiasFrequenciesFile} --fragmentLength {params.fragmentLength} --sampleSize {params.sampleSize} {params.blackListFileName} {input.blackListFileName} --numberOfProcessors {threads} --biasPlot {output.biasPlot} --regionSize {params.regionSize} 2> {log}"

rule correctGCBias:
    params:
        effectiveGenomeSize = get_effectiveGenomeSize,
        binSize = config["correctGCBias"]["binSize"]
    input:
        bamfile = "results/{sample}.filtered.sorted.bam",
        bamfile_index = "results/{sample}.filtered.sorted.bam.bai",
        genome = get_genome(),
        GCbiasFrequenciesFile = "results/{sample}.GCfreq.txt"
    output:
        correctedFile = "results/{sample}.corrected.filtered.sorted.bam"
    log: "results/log/correctGCBias{sample}.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "correctGCBias --bamfile {input.bamfile} --effectiveGenomeSize {params.effectiveGenomeSize} --genome {input.genome} --GCbiasFrequenciesFile {input.GCbiasFrequenciesFile} --correctedFile {output.correctedFile} --binSize {params.binSize} --numberOfProcessors {threads} 2> {log}"

rule multiBamSummary:
    params:
        samples = get_all_samples,
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else ""
    input:
        bamfiles = lambda wildcards: get_all_sorted_bam(wildcards, suffix=".corrected.filtered.sorted.bam", is_results=True),
        bamfiles_index = lambda wildcards: get_all_sorted_bam(wildcards, suffix=".corrected.filtered.sorted.bam.bai", is_results=True),
        BED = get_summary_bed,
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    output:
        outFileName = "results/multiBamSummary.{bed}.npz",
        outRawCounts = "results/multiBamSummary.{bed}.tsv",
        scalingFactors = "results/multiBamSummary.{bed}.scalingFactors"
    log: "results/log/multiBamSummary.{bed}log"
    conda: "../envs/deeptools.yaml"
    shell:
        "multiBamSummary BED-file --bamfiles {input.bamfiles} --outFileName {output.outFileName} --BED {input.BED} --labels {params.samples} {params.blackListFileName} {input.blackListFileName} --numberOfProcessors {threads} --outRawCounts {output.outRawCounts} --scalingFactors {output.scalingFactors} 2> {log}"

rule plotCorrelation:
    params:
        corMethod = config["plotCorrelation"]["corMethod"],
        whatToPlot = config["plotCorrelation"]["whatToPlot"],
        skipZeros = config["plotCorrelation"]["skipZeros"],
        removeOutliers = config["plotCorrelation"]["removeOutliers"]
    input:
        corData = "results/multiBamSummary.{bed}.npz"
    output:
        plotFile = "results/plotCorrelation.{bed}.pdf",
        outFileCorMatrix = "results/plotCorrelation.{bed}.tsv"
    log: "results/log/plotCorrelation.{bed}.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "plotCorrelation --corData {input.corData} --corMethod {params.corMethod} --whatToPlot {params.whatToPlot} --plotFile {output.plotFile} {params.skipZeros} {params.removeOutliers} --outFileCorMatrix {output.outFileCorMatrix} --colorMap bwr --plotNumbers 2> {log}"

rule plotPCA:
    params:
        rowCenter = config["plotPCA"]["rowCenter"],
        ntop = config["plotPCA"]["ntop"]
    input:
        corData = "results/multiBamSummary.{bed}.npz"
    output:
        plotFile = "results/plotPCA.{bed}.pdf",
        outFileNameData = "results/plotPCA.{bed}.tsv"
    log: "results/log/plotPCA.{bed}.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "plotPCA {params.rowCenter} --corData {input.corData} --plotFile {output.plotFile} --outFileNameData {output.outFileNameData} --ntop {params.ntop} 2> {log}"

rule bamCoverage:
    params:
        binSize = config["bamCoverage"]["binSize"],
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else "",
        effectiveGenomeSize = get_effectiveGenomeSize,
        normalizeUsing = config["bamCoverage"]["normalizeUsing"],
        exactScaling = config["bamCoverage"]["exactScaling"],
        skipNonCoveredRegions = config["bamCoverage"]["skipNonCoveredRegions"]
    input:
        bam = "results/{sample}.corrected.filtered.sorted.bam",
        bam_index = "results/{sample}.corrected.filtered.sorted.bam.bai",
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    output:
        outFileName = r"results/{sample,.+(?<!\.averaged)$}.bigwig"
    log: "results/log/bamCoverage.{sample}.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "bamCoverage --bam {input.bam} --outFileName {output.outFileName} --outFileFormat bigwig --binSize {params.binSize} {params.blackListFileName} {input.blackListFileName} --numberOfProcessors {threads} --effectiveGenomeSize {params.effectiveGenomeSize} --normalizeUsing {params.normalizeUsing} {params.exactScaling} {params.skipNonCoveredRegions} 2> {log}"

rule bigwigAverage:
    params:
        binSize = config["bigwigAverage"]["binSize"],
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else ""
    input:
        bigwigs = lambda wildcards: get_group_sorted_bam(wildcards, suffix=".bigwig"),
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    output:
        outFileName = "results/{group}.averaged.bigwig"
    log: "results/log/bigwigAverage.{group}.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "bigwigAverage --bigwigs {input.bigwigs} --binSize {params.binSize} {params.blackListFileName} {input.blackListFileName} --numberOfProcessors {threads} --outFileName {output.outFileName} --outFileFormat bigwig 2> {log}"

rule bigwigCompare:
    params:
        pseudocount = config["bigwigCompare"]["pseudocount"],
        skipZeroOverZero = config["bigwigCompare"]["skipZeroOverZero"],
        operation = config["bigwigCompare"]["operation"],
        skipNonCoveredRegions = config["bigwigCompare"]["skipNonCoveredRegions"],
        binSize = config["bigwigCompare"]["binSize"],
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else ""
    input:
        bigwig1 = "results/{group}.averaged.bigwig",
        bigwig2 = "results/{control_group}.averaged.bigwig",
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    output:
        outFileName = "results/{group}_cmp_{control_group}.bigwig",
    log: "results/log/bigwigCompare.{group}_cmp_{control_group}.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "bigwigCompare --bigwig1 {input.bigwig1} --bigwig2 {input.bigwig2} --pseudocount {params.pseudocount} {params.skipZeroOverZero} --operation {params.operation} {params.skipNonCoveredRegions} --binSize {params.binSize} {params.blackListFileName} {input.blackListFileName} --numberOfProcessors {threads} --outFileName {output.outFileName} --outFileFormat bigwig 2> {log}"

rule computeMatrix:
    params:
        referencePoint = config["computeMatrix"]["referencePoint"],
        beforeRegionStartLength = config["computeMatrix"]["beforeRegionStartLength"],
        afterRegionStartLength = config["computeMatrix"]["afterRegionStartLength"],
        nanAfterEnd = config["computeMatrix"]["nanAfterEnd"],
        binSize = config["computeMatrix"]["binSize"],
        missingDataAsZero = config["computeMatrix"]["missingDataAsZero"],
        skipZeros = config["computeMatrix"]["skipZeros"],
        minThreshold = config["computeMatrix"]["minThreshold"],
        maxThreshold = config["computeMatrix"]["maxThreshold"],
        blackListFileName = "--blackListFileName" if config["rmsk_url"] else ""
    input:
        regionsFileName = config["computeMatrix"]["regionsFileName"],
        scoreFileName = get_all_bigwig,
        blackListFileName = "results/rmsk.bed" if config["rmsk_url"] else []
    output:
        outFileName = "results/computeMatrix.gz",
        outFileNameMatrix = "results/computeMatrix.tsv",
    log: "results/log/computeMatrix.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "computeMatrix reference-point --regionsFileName {input.regionsFileName} --scoreFileName {input.scoreFileName} --outFileName {output.outFileName} --outFileNameMatrix {output.outFileNameMatrix} --referencePoint {params.referencePoint} --beforeRegionStartLength {params.beforeRegionStartLength} --afterRegionStartLength {params.afterRegionStartLength} {params.nanAfterEnd} --binSize {params.binSize} {params.missingDataAsZero} {params.skipZeros} {params.minThreshold} {params.maxThreshold} {params.blackListFileName} {input.blackListFileName} --smartLabels --numberOfProcessors {threads} 2> {log}"

rule plotHeatmap:
    params:
        interpolationMethod = config["plotHeatmap"]["interpolationMethod"],
        dpi = config["plotHeatmap"]["dpi"],
        plotType = config["plotHeatmap"]["plotType"],
        sortRegions = config["plotHeatmap"]["sortRegions"],
        sortUsing = config["plotHeatmap"]["sortUsing"],
        sortUsingSamples = config["plotHeatmap"]["sortUsingSamples"],
        averageTypeSummaryPlot = config["plotHeatmap"]["averageTypeSummaryPlot"],
        zMin = config["plotHeatmap"]["zMin"],
        zMax = config["plotHeatmap"]["zMax"],
        heatmapHeight = config["plotHeatmap"]["heatmapHeight"],
        heatmapWidth = config["plotHeatmap"]["heatmapWidth"],
        xAxisLabel = config["plotHeatmap"]["xAxisLabel"],
        yAxisLabel = config["plotHeatmap"]["yAxisLabel"],
        refPointLabel = config["plotHeatmap"]["refPointLabel"],
        regionsLabel = get_all_regionsFileName_stem,
        samplesLabel = get_all_bigwig_stem,
        yMin = config["plotHeatmap"]["yMin"],
        yMax = config["plotHeatmap"]["yMax"],
        legendLocation = config["plotHeatmap"]["legendLocation"],
        perGroup = config["plotHeatmap"]["perGroup"]
    input:
        matrixFile = "results/computeMatrix.gz"
    output:
        outFileName = "results/plotHeatmap.pdf",
        outFileSortedRegions = "results/plotHeatmap.bed"
    log: "results/log/plotHeatmap.log"
    conda: "../envs/deeptools.yaml"
    shell:
        "plotHeatmap --matrixFile {input.matrixFile} --outFileName {output.outFileName} --outFileSortedRegions {output.outFileSortedRegions} --interpolationMethod {params.interpolationMethod} --dpi {params.dpi} --plotType {params.plotType} --sortRegions {params.sortRegions} --sortUsing {params.sortUsing} --sortUsingSamples {params.sortUsingSamples} --averageTypeSummaryPlot {params.averageTypeSummaryPlot} --zMin {params.zMin} --zMax {params.zMax} --heatmapHeight {params.heatmapHeight} --heatmapWidth {params.heatmapWidth} --xAxisLabel {params.xAxisLabel} --yAxisLabel {params.yAxisLabel} --refPointLabel {params.refPointLabel} --regionsLabel {params.regionsLabel} --samplesLabel {params.samplesLabel} {params.yMin} {params.yMax} --legendLocation {params.legendLocation} --perGroup {params.perGroup} 2> {log}"


