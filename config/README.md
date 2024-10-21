# Instruction
effectiveGenomeSize_normsk.csv and effectiveGenomeSize_rmsk.csv are copied from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html. They are used to esitimate the effective genome size after remove repeats or not.

config.yaml contains settings for deeptools commands. Their meanings are described in https://deeptools.readthedocs.io/en/latest/.

samplesheet.csv contains four columns:
```list
sample is the file name of sorted bam.
group marks samples as replicates, say replicates of WT.
control_group indicates the control of the current group.
peak is the peaks called by MACS(2|3) from the sample bam.
```