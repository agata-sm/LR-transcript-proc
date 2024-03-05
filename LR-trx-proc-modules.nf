// modules for minimal.nf


params.st="stringtie"
params.stOut="${params.outdir}/${params.st}"

params.stm="stringtie_merge"
params.stmOut="${params.outdir}/${params.stm}"

params.stm_gffcmp="stringtie_gff_compare"
params.stGffCmpOut="${params.outdir}/${params.stm_gffcmp}"

params.esp="ESPRESSO"
params.espOut="${params.outdir}/${params.esp}"

params.espm="espresso_merge"
params.espmOut="${params.outdir}/${params.espm}"

params.espm_gffcmp="espresso_gff_compare"
params.espGffCmpOut="${params.outdir}/${params.espm_gffcmp}"




// assets
//params.countertemplate="${projectDir}/assets/template.properties"


// scripts
params.scripts="${projectDir}/bin"

// versions
params.verfile="software.versions"

/// modules


process stringtie {
    publishDir params.stOut, mode:'copy'

    label 'mid_mem'

    input:
    tuple path(bamfile), val(smpl_id)

    output:
    path("${smpl_id}.stringtie.gtf"), emit: stringtie_gtf_ch

    script:
    """
    echo $bamfile
    echo $smpl_id
    
    stringtie $bamfile --rf -L -p ${task.cpus} -v -G ${params.refGTF} >${smpl_id}.stringtie.gtf

    date >>${params.verfile}
    echo "stringtie" >>${params.verfile}
    stringtie --version >>${params.verfile}
    echo "" >>${params.verfile}
    """

}

process stringtie_merge {
    publishDir params.stmOut, mode:'copy'

    label 'mid_mem'

    input:
    path stringtie_gtfs

    output:
    path("${params.projname}.stringtie.merged.gtf"), emit: stringtie_merged_ch

    script:
    """
    echo ${stringtie_gtfs}

    stringtie --merge ${stringtie_gtfs}  -G ${params.refGTF} -o ${params.projname}.stringtie.merged.gtf

    date >>${params.verfile}
    echo "stringtie" >>${params.verfile}
    stringtie --version >>${params.verfile}
    echo "" >>${params.verfile}
    """

}


process gffcompare_stringtie {
    publishDir params.stGffCmpOut, mode:'copy'

    label 'small'

    input:
    path stringtie_merged

    output:
    path("gffcompare_stringtie.stats")
    path("gffcompare_stringtie.annotated.gtf")
    path("gffcompare_stringtie.loci")
    path("gffcompare_stringtie.loci.${params.projname}.stringtie.merged.gtf.refmap")
    path("gffcompare_stringtie.loci.${params.projname}.stringtie.merged.gtf.tmap")


    script:
    """
    echo ${stringtie_merged}

    gffcompare -R -r ${params.refGTF} -o gffcompare_stringtie $stringtie_merged
    
    date >>${params.verfile}
    echo "gffcompare" >>${params.verfile}
    gffcompare --version >>${params.verfile}
    echo "" >>${params.verfile}
    """

}



process espresso {
    publishDir params.espOut, mode:'copy'

    label 'big_mem'

    input:
    tuple path(bamfile), val(smpl_id)
    
    output:
    path "${smpl_id}/espresso_s_summary.txt"
    path "${smpl_id}/SJ_group_all.fa"
    path "${smpl_id}/*abundance.esp"
    path "${smpl_id}/*updated.gtf", emit: espresso_gtf_ch
    path "${smpl_id}/1_SJ_simplified.list"
    path "${smpl_id}/espresso_q_summary.txt"
    path "${smpl_id}/0/1_read_final.txt"
    path "${smpl_id}/0/espresso_c_summary.txt"
    path "${smpl_id}/0/sam.list3"
    path "${smpl_id}/0/sj.list"


    script:
    """
    mkdir ${smpl_id}

    echo "${bamfile}\t${smpl_id}" >>sample_sheet.tsv

    #perl /espresso/src/ESPRESSO_S.pl -M M -L sample_sheet.tsv -T ${task.cpus} --sort_buffer_size ${task.memory} -F ${params.refFa} -A ${params.refGTF} -O ${smpl_id}
    #perl /espresso/src/ESPRESSO_C.pl -I ${smpl_id} -F ${params.refFa} -X 0 -T ${task.cpus} --sort_buffer_size ${task.memory}
    #perl /espresso/src/ESPRESSO_Q.pl -L "${smpl_id}/sample_sheet.tsv.updated" -A ${params.refGTF} -T ${task.cpus}

    perl /espresso/src/ESPRESSO_S.pl -M M -L sample_sheet.tsv -T ${params.espresso_cpus} --sort_buffer_size ${params.espresso_mem} -F ${params.refFa} -A ${params.refGTF} -O ${smpl_id}
    perl /espresso/src/ESPRESSO_C.pl -I ${smpl_id} -F ${params.refFa} -X 0 -T ${params.espresso_cpus} --sort_buffer_size ${params.espresso_mem}
    perl /espresso/src/ESPRESSO_Q.pl -L "${smpl_id}/sample_sheet.tsv.updated" -A ${params.refGTF} -T ${params.espresso_cpus}


    date >>${params.verfile}
    echo "ESPRESSO" >>${params.verfile}
    perl /espresso/src/ESPRESSO_S.pl --help >>${params.verfile}
    echo "" >>${params.verfile}
    perl /espresso/src/ESPRESSO_C.pl --help >>${params.verfile}
    echo "" >>${params.verfile}
    perl /espresso/src/ESPRESSO_Q.pl --help >>${params.verfile}
    echo "" >>${params.verfile}
    """
}


process espresso_merge {
    publishDir params.espmOut, mode:'copy'

    label 'mid_mem'

    input:
    path espresso_gtfs

    output:
    path("${params.projname}.espresso.merged.gtf"), emit: espresso_merged_ch

    script:
    """
    echo ${espresso_gtfs}

    stringtie --merge ${espresso_gtfs}  -G ${params.refGTF} -o ${params.projname}.espresso.merged.gtf

    date >>${params.verfile}
    echo "stringtie" >>${params.verfile}
    stringtie --version >>${params.verfile}
    echo "" >>${params.verfile}
    """

}


process gffcompare_espresso {
    publishDir params.espGffCmpOut, mode:'copy'

    label 'small'

    input:
    path espresso_merged

    output:
    path("gffcompare_espresso.stats")
    path("gffcompare_espresso.annotated.gtf")
    path("gffcompare_espresso.loci")
    path("gffcompare_espresso.loci.${params.projname}.espresso.merged.gtf.refmap")
    path("gffcompare_espresso.loci.${params.projname}.espresso.merged.gtf.tmap")


    script:
    """
    echo ${espresso_merged}

    gffcompare -R -r ${params.refGTF} -o gffcompare_espresso $espresso_merged
    
    date >>${params.verfile}
    echo "gffcompare" >>${params.verfile}
    gffcompare --version >>${params.verfile}
    echo "" >>${params.verfile}
    """

}



