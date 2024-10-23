process stringtie {
    label 'stringtie'
    label 'mid_mem'

    tag {smpl_id}

    container = 'docker://ctomlins/stringtie2'

    input:
    tuple path(bamfile), val(smpl_id)

    output:
    path("*.gtf"), emit: stringtie_gtf_ch
    path "versions.txt"

    script:

    def args = task.ext.args ?: ''

    """
    stringtie $bamfile \\
        ${args} \\
        -p ${task.cpus} \\
        -G ${params.refGTF} >${smpl_id}.stringtie.gtf

    cat <<-END_VERSIONS > versions.txt
    Software versions for LR-trx-proc.nf
    \$( date )
    process ** stringtie **
    stringtie
    \$( stringtie --version )
    END_VERSIONS
    """

}