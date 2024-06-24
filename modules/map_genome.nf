process map_genome {

    label 'mid_mem'
    label 'wftrx'
    container = 'docker://nanozoo/minimap2:2.28--9e3bd01'

    tag {smpl_id}

    input:
    tuple path(fastq_proc), val(smpl_id), path(genome_idx_ch)

    output:
    tuple path("${smpl_id}_all_alns.minimap2.bam"), val(smpl_id), emit: mapped_genome_ch
    //tuple path("${smpl_id}_all_alns.minimap2.bam"), path("${smpl_id}_all_alns.minimap2.bam.bai"), val(smpl_id), emit: mapped_genome_ch
    path "versions.txt"

    script:

    def args = task.ext.args ?: ''

    """
    minimap2 -t ${task.cpus} \\
    ${args} -a \\ 
    ${genome_idx_ch} ${fastq_proc} \\
    | samtools sort -@ ${task.cpus} \\
    --write-index \\ 
    -o "${smpl_id}_all_alns.minimap2.bam" -

    cat <<-END_VERSIONS > versions.txt
    Software versions for LR-trx-proc.nf
    \$( date )
    process ** map_genome **
    minimap2
    \$( minimap2 --version )
    samtools
    \$( samtools --version )
    END_VERSIONS
    """


}