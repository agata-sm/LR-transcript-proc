#! /usr/bin/env nextflow


/* 
 * pipeline for analysis of ONT LR transcriptomics data (cDNA-seq)
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
params.projdir = "${launchDir}/${params.projname}"
params.outdir = "${params.projdir}/${params.resdir}"


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

include { preprocess_reads; genome_idx; map_genome; stringtie; stringtie_merge; gffcompare_stringtie; gffcompare_espresso; espresso_s_input; espresso_c_smpl; espresso_q_input} from './LR-trx-proc-modules.nf'


///////////////////////////
//channels

println ""

println "Sample paths"


//samples channel
smpls_ch= Channel.fromPath(params.samplesheet, checkIfExists:true)
	smpls_ch
	    .splitCsv(header:true, sep: '\t', strip: true)
	    .map{ row-> tuple(row.path, row.sample) }
	    .view()
	    .set { smpls_ch }

println ""

println "Genome reference"

//genome fa channel
genome_ch= Channel.fromPath(params.refFa, checkIfExists:true)
	.view()

println ""




/////////////////////////////
// workflows

//default

workflow {

	//cat fastq files and preprocess using pychopper
	preprocess_reads(smpls_ch)

	//minimap2
	genome_idx(genome_ch)

	full_len_reads_map_ch=preprocess_reads.out.full_len_reads
	full_len_reads_map_ch
			.combine(genome_idx.out.genome_idx_ch)
			.set {full_len_reads_map_ch}
	map_genome(full_len_reads_map_ch)

	//stringtie
	stringtie(map_genome.out.mapped_genome_ch)
	stringtie_out_ch=stringtie.out.stringtie_gtf_ch
	stringtie_merge(stringtie_out_ch.collect())
	gffcompare_stringtie(stringtie_merge.out.stringtie_merged_ch)

	//espresso the right way
	mapped_genome_all_ch=map_genome.out.mapped_genome_ch
		.collect()

	// ESPRESSO
	espresso_s_input(mapped_genome_all_ch)
	espresso_s_out_ch=espresso_s_input.out.espresso_s_out_ch

	// map sample ID to espresso S index
	sample_idx_ch= espresso_s_input.out.espresso_s_samplesheet_ch
	sample_idx_ch
		    .splitCsv(header:false, sep: '\t', strip: true)
		    .map{ row -> tuple( "${row[1]}","${row[2]}") }
		    .set { sample_idx_ch }

	espresso_c_smpl(sample_idx_ch, espresso_s_out_ch)
	espresso_c_for_q_ch=espresso_c_smpl.out.espresso_c_smpl_ch
	espresso_q_input(espresso_c_for_q_ch.collect(),espresso_s_input.out.espresso_s_samplesheet_ch)

	gffcompare_espresso(espresso_q_input.out.espresso_gtf_ch)


	// SQANTI



}



