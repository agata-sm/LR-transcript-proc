process espresso_q_input {

    label 'espresso_q_input'

    container = 'docker://agatasm/espresso-1.4.0-perl'

    input:
    path espresso_c_smpl_ch
    path espresso_s_samplesheet_ch

    output:
    path "*updated.gtf", emit: espresso_gtf_ch
    path "*abundance.esp"
    path "versions.txt"

    script:
    def ext_args = task.ext.args ?: ''
    """
    perl /espresso/src/ESPRESSO_Q.pl  ${ext_args} \\
        -L sample_sheet.tsv.updated \\
        -A ${params.refGTF} \\
        -T ${cpus} \\
        -O .

    cat <<-END_VERSIONS > versions.txt
    Software versions for LR-trx-proc.nf
    \$( date )
    process ** espresso_q_input **
    ESPRESSO_Q.pl
    \$(perl /espresso/src/ESPRESSO_Q.pl --help | grep -E "Program | Version" )
    END_VERSIONS
    """


}