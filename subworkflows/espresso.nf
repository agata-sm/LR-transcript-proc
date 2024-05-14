/*
subworklow to process LR mapping to genome for isoform detection and quantification
takes bam files and runs ESPRESSO (DOI: 10.1126/sciadv.abq5072, https://github.com/Xinglab/espresso)
*/


params.espressoOut="${params.outdir}/espresso"
params.scriptDir="${projectDir}/bin"

params.espressoOutS="${params.espressoOut}/S"


process espresso_s {

    publishDir params.espressoOutS, mode:'copy',
    saveAs: {filename ->
        if (filename.endsWith(".gtf")) "$filename"
        else if (filename.endsWith(".tsv")) "$filename"
        else if (filename.endsWith(".txt")) "$filename"
        else null
    }

    label 'espressoS'

    input:
    tuple path(bamfile), val(smpl_id)

    output:
    path("*.tsv")
    path "versions.txt"

    script:
    """
    echo "${bamfile}\t${smpl_id}" >>sample_sheet.tsv


    #perl /espresso/src/ESPRESSO_S.pl -M M -L sample_sheet.tsv -T ${params.espresso_cpus} --sort_buffer_size ${params.espresso_mem} -F ${params.refFa} -A ${params.refGTF} -O ${smpl_id}


    cat <<-END_VERSIONS > versions.txt
    Software versions for LR-trx-proc.nf
    \$( date )
    process ** espresso_s **
    ESPRESSO
    \$(  perl /espresso/src/ESPRESSO_S.pl --help | head -n 5 )
    END_VERSIONS
    """

}


workflow espresso {
    take:
       smpls_ch

    main:
    	espresso_s(smpls_ch)

}


