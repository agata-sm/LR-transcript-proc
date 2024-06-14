# LR-transcript-proc
nextflow pipeline for analysis of ONT cDNA sequencing data

# Test run on Rackham

at `/sw/courses/epigenomics/agata/lr/runs/LR_pipeline_14vi/`

cmds:

```
    module load bioinfo-tools
    module load Nextflow/22.10.1

    export NXF_HOME="/sw/courses/epigenomics/agata/lr/software/nf"
    export APPTAINERENV_TMPDIR="/proj/naiss2023-23-205/nobackup/private/nbis6556/tsts/LR-pipeline/containers"
    export NXF_SINGULARITY_CACHEDIR="/proj/naiss2023-23-205/nobackup/private/nbis6556/tsts/LR-pipeline/containers"

    export pipelineDir="/sw/courses/epigenomics/agata/lr/runs/LR_pipeline_14vi/gh/LR-transcript-proc"


		nextflow run ${pipelineDir}/LR-trx-proc.nf -c proj-LR-tst.config -profile singularity,cluster
```

containers at `/proj/naiss2023-23-205/nobackup/private/nbis6556/tsts/LR-pipeline/containers` 
and `/proj/naiss2023-23-205/nobackup/private/nbis6556/software/containers/` 
and `/sw/courses/epigenomics/agata/lr/software/containers/`


