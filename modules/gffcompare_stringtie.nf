process gffcompare_stringtie {

    label 'small'
    container = 'docker://agatasm/gffcompare-0.12.9'

    input:
    path stringtie_merged

    output:
    path("gffcompare_stringtie.stats")
    path("gffcompare_stringtie.annotated.gtf")
    path("gffcompare_stringtie.loci")
    path("gffcompare_stringtie.${stringtie_merged}.refmap")
    path("gffcompare_stringtie.${stringtie_merged}.tmap")
    path "versions.txt"

    script:
    def args = task.ext.args ?: ''
    """
    gffcompare ${args} \\
        -r ${params.refGTF} \\
        -o gffcompare_stringtie $stringtie_merged
    
    cat <<-END_VERSIONS > versions.txt
    Software versions for LR-trx-proc.nf
    \$( date )
    process **  gffcompare_stringtie **
    gffcompare
    \$( gffcompare --version )
    END_VERSIONS
    """

}