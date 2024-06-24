process espresso_c_smpl {
    
    label 'espresso_c_smpl'

    container = 'docker://agatasm/espresso-1.4.0:perl'

    tag {smpl_id}

    input:
    tuple val(smpl_id), val(smpl_idx)
    path espresso_s_out_ch

    output:
    path "espressoS/${smpl_idx}", emit: espresso_c_smpl_ch
    path "${smpl_id}.espresso_c_summary.txt"
    path "versions.txt"


    script:
    def ext_args = task.ext.args ?: ''
    """
    cp -r ESPRESSO_S espressoS ## to isolate the execution process

    perl /espresso/src/ESPRESSO_C.pl ${ext_args} \\
        -I espressoS \\
        -T ${task.cpus} \\
        --sort_buffer_size ${task.memory} \\
        -F ${params.refFa} \\
        -X ${smpl_idx}

    cp "espressoS/${smpl_idx}/espresso_c_summary.txt" "${smpl_id}.espresso_c_summary.txt"

    cat <<-END_VERSIONS > versions.txt
    Software versions for LR-trx-proc.nf
    \$( date )
    process **  espresso_c_smpl **
    ESPRESSO_C.pl
    \$(perl /espresso/src/ESPRESSO_C.pl --help | grep -E "Program | Version" )
    END_VERSIONS

    """

}