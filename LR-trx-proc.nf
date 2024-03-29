#! /usr/bin/env nextflow


/* 
 * pipeline for analysis of LR transcriptomics data processed by wf-transcriptomes (tested with our adapted fork)
 * Nextflow DSL2
 * 
 * Author: Agata Smialowska
 * January 2024
 */ 


nextflow.enable.dsl=2

params.pipelinename="LR transcriptomics analysis pipeline"

/* 
 * pipeline input parameters 
 */
params.resdir = "results"
params.projdir = "$launchDir/${params.projname}"

params.outdir = "${params.projdir}/${params.resdir}"

params.logdir = 'logs'


log.info """\
 LONG READ TRANSCRIPTOME ANALYSIS - N F   P I P E L I N E
 ===================================
 
 sample sheet: ${params.samplesheet}

 outdir       : ${params.outdir}
 """
 .stripIndent()

println ""

println ""
println ""

///////////////////////////
//channels

//samples channel
smpls_ch= Channel.fromPath(params.samplesheet, checkIfExists:true)
	smpls_ch
	    .splitCsv(header:true, sep: '\t', strip: true)
	    .map{ row-> tuple(row.path, row.sample) }
	    //.view()
	    .set { smpls_ch }

println ""



/////////////////////////////
// processes
include { stringtie; espresso; stringtie_merge; gffcompare_stringtie; espresso_merge; gffcompare_espresso} from './LR-trx-proc-modules.nf'



/////////////////////////////
// workflows

//default

workflow {

	//stringtie
	stringtie(smpls_ch)
	stringtie_out_ch=stringtie.out.stringtie_gtf_ch
	stringtie_merge(stringtie_out_ch.collect())
	gffcompare_stringtie(stringtie_merge.out.stringtie_merged_ch)

	//espresso
	espresso(smpls_ch)
	espresso_out_ch=espresso.out.espresso_gtf_ch
	espresso_merge(espresso_out_ch.collect())
	gffcompare_espresso(espresso_merge.out.espresso_merged_ch)

}


