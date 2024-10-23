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

### Some feedback for the 1st run on the dardel compute cluster

The average compute node has 256 cores but only 232Gb RAM. As the memory is evenly distributed
among the 256 cores, this results in only 0.9Gb per core! And the cores are updated according to the requested memory for every slurm jobs, which usually results in a high penalty of "unused CPUs".
Therefore we need to get a good estimate of the required memory for each step of the pipeline.
And if a memory increase for a process is required, then try to increase the CPU-cores as well (as you get charged anyway for those).

**0. Data**

- human reference genome
- cDNA data set:

| barcode |    num_seq |        sum_len |
| :------ | ---------: | -------------: |
| 01      | 22,566,774 | 21,752,143,183 |
| 02      | 23,179,382 | 19,244,463,129 |

**1. recommended setup for next runs with human reference genome on Dardel compute cluster:**

```bash
    withName: 'preprocess_reads' {
        memory = 40.GB
        cpus   = 40
        time   = 6.h
    }
    withName: 'genome_idx' {
        memory = 18.GB
        cpus   = 20
        time   = 2.h
    }
    withName: 'map_genome' {
        memory = 30.GB
        cpus   = 30
        time   = 6.h
    }
    withName: 'stringtie' {
        memory = 6.GB
        cpus   = 8
        time   = 2.h
    }
    withName: 'stringtie_merge' {
        memory = 1.GB
        cpus   = 1
        time   = 1.h
    }
    withName: 'gffcompare_stringtie' {
        memory = 2.GB
        cpus   = 1
        time   = 1.h
    }
    withName: 'espresso_s_input' {
        memory = 30.GB
        cpus   = 30
        time   = 6.h

    }
    withName: 'espresso_c_smpl' {
        memory = 45.GB
        cpus   = 45
        time   = 2.d
    }
    withName: 'espresso_q_input' {
        memory = 140.GB
        cpus   = 80
        time   = 6.h
    }
    withName: 'gffcompare_espresso' {
        memory = 2.GB
        cpus   = 1
        time   = 1.h
    }
    withName: 'sqanti_qc' {
        memory = 8.GB
        cpus   = 10
        time   = 2.h
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
