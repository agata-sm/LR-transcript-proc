# LR-transcript-proc

nextflow pipeline for analysis of ONT cDNA sequencing data

# Update

- separate modules into individual files
- add uppmax and dardel config files
- all program arguments can be set in modules.config
- the pipeline uses the following docker container:
  - [docker://agatasm/espresso-1.4.0:perl](docker://agatasm/espresso-1.4.0:perl)
  - [docker://nanozoo/minimap2:2.28--9e3bd01](docker://nanozoo/minimap2:2.28--9e3bd01) (NEW)
  - [docker://agatasm/gffcompare:0.12.9](docker://agatasm/gffcompare:0.12.9)
  - [docker://ontresearch/wf-transcriptomes:sha203915eb4b4dd444cb2e845d0b9f7814e26b7b5c](docker://ontresearch/wf-transcriptomes:sha203915eb4b4dd444cb2e845d0b9f7814e26b7b5c)
  - [docker://pippel/squanti3:v5.2.1](docker://pippel/squanti3:v5.2.1) (NEW)
  - [docker://ctomlins/stringtie2](docker://ctomlins/stringtie2)

The pipeline can be run on uppmax by specifying the uppmax profile:

```bash
# Set common path to store all Singularity containers
STORAGEALLOC="/proj/snic2022-6-208"
export NXF_SINGULARITY_CACHEDIR="${STORAGEALLOC}/nobackup/ebp-singularity-cache"
export pipelineDir="/proj/snic2022-6-208/private/martin/support_7486/gh/LR-transcript-proc"

set -euo pipefail

# Activate shared Nextflow environment
eval "$(conda shell.bash hook)"
conda activate "${STORAGEALLOC}/conda/nextflow-env"

export NXF_VER=23.10.1

nextflow run ${pipelineDir}/LR-trx-proc.nf \
        -c nextflow.config \
        -profile uppmax \
        -params-file workflow_parameters.yml \
        -ansi-log false \
        -resume
```

... or on dardel by specifying the dardel profile:

```bash
# Set common path to store all Singularity containers
STORAGEALLOC="/cfs/klemming/projects/snic/naiss2024-23-311/martin/"
export NXF_SINGULARITY_CACHEDIR="${STORAGEALLOC}/cache"
export pipelineDir="/cfs/klemming/projects/snic/naiss2024-23-311/martin/LR_pipeline/gh/LR-transcript-proc-myFORK"

set -euo pipefail

# Activate shared Nextflow environment
eval "$(/sw/apps/conda/latest/rackham_stage/condabin/conda shell.bash hook)"
conda activate /cfs/klemming/scratch/p/pippel/prog/conda_envs/nextflow

export NXF_VER=23.10.1

nextflow run ${pipelineDir}/LR-trx-proc.nf \
        -c nextflow.config \
        -profile dardel \
        -params-file workflow_parameters.yml \
        -ansi-log false \
        -resume
```

The `nextflow.config` file looks like this:

```bash
params {
        projname      = "KLK4-cdna_bc1-bc2"
        project       = "naiss2024-22-675"

        refFa         = "/cfs/klemming/projects/snic/naiss2024-23-311/martin/ref/gencode_GRCh38.p14/GRCh38.primary_assembly.genome.fa"
        refGTF        = "/cfs/klemming/projects/snic/naiss2024-23-311/martin/ref/gencode_GRCh38.p14/gencode.v46.primary_assembly.annotation.gtf"

        samplesheet   = "sample_paths.tsv"

        espresso_mito = "chrM"
}
```

# Original README:

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
