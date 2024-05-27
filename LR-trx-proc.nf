#! /usr/bin/env nextflow


/* 
 * pipeline for analysis of LR transcriptomics data processed by wf-transcriptomes (tested with our adapted fork)
 * Nextflow DSL2
 * 
 * Author: Agata Smialowska
 * January - May 2024
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


/////////////////////////////
// processes
//include { stringtie; espresso_run; stringtie_merge; gffcompare_stringtie; espresso_merge; gffcompare_espresso; map_genome} from './LR-trx-proc-modules.nf'

include { map_genome; stringtie; espresso_run; stringtie_merge; gffcompare_stringtie; espresso_merge; gffcompare_espresso; espresso_s_input; espresso_c_smpl; espresso_q_input} from './LR-trx-proc-modules.nf'


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
// workflows

//default

workflow {

	//stringtie
	stringtie(smpls_ch)
	stringtie_out_ch=stringtie.out.stringtie_gtf_ch
	stringtie_merge(stringtie_out_ch.collect())
	gffcompare_stringtie(stringtie_merge.out.stringtie_merged_ch)

	//espresso
	//espresso(smpls_ch)
	//espresso_out_ch=espresso.out.espresso_gtf_ch
	//espresso_merge(espresso_out_ch.collect())
	//gffcompare_espresso(espresso_merge.out.espresso_merged_ch)

	//mock process atm but later it will actually map reads to genome
	map_genome(smpls_ch)

	//espresso the right way
	mapped_genome_ch=map_genome.out.mapped_genome_ch
		.collect()

	// ESPRESSO
	espresso_s_input(mapped_genome_ch)
	espresso_s_out_ch=espresso_s_input.out.espresso_s_out_ch

	// map sample ID to espresso S index
	sample_idx_ch= espresso_s_input.out.espresso_s_samplesheet_ch
	sample_idx_ch
		    .splitCsv(header:false, sep: '\t', strip: true)
		    //.view( row-> "${row[1]}, ${row[2]}" )
		    .map{ row -> tuple( "${row[1]}","${row[2]}") }
		    .set { sample_idx_ch }

	espresso_c_smpl(sample_idx_ch, espresso_s_out_ch)
	espresso_c_for_q_ch=espresso_c_smpl.out.espresso_c_smpl_ch
	espresso_q_input(espresso_c_for_q_ch.collect(),espresso_s_input.out.espresso_s_samplesheet_ch)

}



