process stringtie_merge {
    label 'stringtie_merge'
    label 'mid_mem'
    container = 'docker://ctomlins/stringtie2'

    input:
    path stringtie_gtfs

    output:
    path("${params.projname}.stringtie.merged.gtf"), emit: stringtie_merged_ch
    path "versions.txt"

    script:
    
    def args = task.ext.args ?: ''
    
    """
    stringtie \\
        ${args} \\
        --merge ${stringtie_gtfs} \\
        -G ${params.refGTF} \\
        -o ${params.prefixOut}.stringtie.merged.gtf

    cat <<-END_VERSIONS > versions.txt
    Software versions for LR-trx-proc.nf
    \$( date )
    process ** stringtie_merge **
    stringtie
    \$( stringtie --version )
    END_VERSIONS

    """

}