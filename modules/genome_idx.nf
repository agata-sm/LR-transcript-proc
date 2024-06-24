process genome_idx {

    label 'mid_mem'
    
    container = 'docker://nanozoo/minimap2:2.28--9e3bd01'

    input:
    path genome_ch

    output:
    path "genome.minimap.idx.mmi", emit: genome_idx_ch

    script:
    
    def args = task.ext.args ?: ''

    """
    minimap2 -t ${task.cpus} ${args} -d "genome.minimap.idx.mmi" ${genome_ch}
    """
}