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
params.resdir    = "results"
params.projdir   = "${launchDir}/${params.projname}"
params.outdir    = "${params.projdir}/${params.resdir}"
params.prefixOut = "${params.projname}"
params.seqtype   = "cDNA"      // possible values ["cDNA", "dRNA"]

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
include { preprocess_reads     } from "$projectDir/modules/preprocess_reads.nf"
include { genome_idx           } from "$projectDir/modules/genome_idx.nf"
include { map_genome           } from "$projectDir/modules/map_genome.nf"
include { stringtie            } from "$projectDir/modules/stringtie.nf"
include { stringtie_merge      } from "$projectDir/modules/stringtie_merge.nf"
include { gffcompare_stringtie } from "$projectDir/modules/gffcompare_stringtie.nf"
include { gffcompare_espresso  } from "$projectDir/modules/gffcompare_espresso.nf"
include { espresso_s_input     } from "$projectDir/modules/espresso_s_input.nf"
include { espresso_c_smpl      } from "$projectDir/modules/espresso_c_smpl.nf"
include { espresso_q_input     } from "$projectDir/modules/espresso_q_input.nf"
include { sqanti_qc            } from "$projectDir/modules/sqanti_qc.nf"



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
	sqanti_qc(espresso_q_input.out.espresso_gtf_ch)

}

workflow.onComplete {
    if( workflow.success ){
        log.info("""
        Thanks you for using the LR-transcript-proc workflow.
        The workflow completed successfully.

        Results are located in the folder: $params.outdir
        """)
    } else {
        log.info("""
        Thanks you for using the LR-transcript-proc workflow.
        The workflow completed unsuccessfully.
        """)
    }
}