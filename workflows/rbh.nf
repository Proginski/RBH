/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WELCOME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// log.info ""
// log.info ">>> RBH workflow <<<"
// log.info ""
// log.info "This workflow takes a multiFASTA (A) and finds its reciprocal best hits among multiFASTAs (B)."
// log.info ""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS MANAGMENT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECK_INPUTS                                      } from '../modules/local/check_inputs.nf'
include { DIAMOND_BLASTP_A_VS_B                             } from '../modules/local/diamond_blastp_A_vs_B.nf'
include { BLAST_BEST_HITS                                   } from '../modules/local/blast_best_hits.nf'
include { FASOMERECORDS_B_BH                                } from '../modules/local/faSomeRecords_B_bh.nf'
include { DIAMOND_BLASTP_A_VS_B as DIAMOND_BLASTP_B_BH_VS_A } from '../modules/local/diamond_blastp_A_vs_B.nf'
include { BLAST_BEST_HITS as BLAST_BEST_HITS_REVERSE        } from '../modules/local/blast_best_hits.nf'
include { RECIPROCAL_BEST_HITS                              } from '../modules/local/reciprocal_best_hits.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RBH {

	// A channel with the B fastas
	if (params.includeA){
		log.info "WARNING '--includeA' is on : searching for the A orthologs might take a while."
		B_faa_ch = Channel.fromPath(params.B)
	} else {
		log.info "Excluding the A fasta (default) (use '--includeA' to force its inclusion)."
		B_faa_ch = Channel.fromPath(params.B)
		.filter { faa -> faa != file(params.A) }
	}
	CHECK_INPUTS(
		file(params.A),
		B_faa_ch.collect(),
		params.fasta_type
	)
	A_faa_ch = CHECK_INPUTS.out.A_faa
	B_faa_ch = CHECK_INPUTS.out.B_faas
		.flatten()

	// BLASTp the query fasta against each of the B fastas
	DIAMOND_BLASTP_A_VS_B(
		A_faa_ch,
		B_faa_ch
	)

	// Get a two col tsv of best hits (query, query_best_hit)
	BLAST_BEST_HITS( DIAMOND_BLASTP_A_VS_B.out )

	// Get a subset of each B fasta with only the entries which were among the best hits
	FASOMERECORDS_B_BH( BLAST_BEST_HITS.out )

	// BLASTp each best_hits fasta against the query fasta
	DIAMOND_BLASTP_B_BH_VS_A( 
		FASOMERECORDS_B_BH.out
		.map { A_faa, B_faa, B_bh_faa -> B_bh_faa },
		A_faa_ch
	)

	// Get the best hits tsv of DIAMOND_BLASTP_B_BH_VS_A
	BLAST_BEST_HITS_REVERSE( DIAMOND_BLASTP_B_BH_VS_A.out )

	bh_ch =  BLAST_BEST_HITS.out
	.map{ A_faa, B_faa, bh -> [ B_faa.getBaseName(), bh ] }

	rev_bh_ch = BLAST_BEST_HITS_REVERSE.out
	.map{ A_faa, B_faa, bh -> [ file(file(A_faa.getBaseName()).getBaseName()).getBaseName(), bh ] } // Here there are three extensions to remove, and getSimpleName() would damage names that contain dots.

	recip_bh_ch = bh_ch.combine(rev_bh_ch, by: 0)
	.map{ B_raw, bh_AB, bh_BA -> [ file(params.A).getBaseName(), B_raw, bh_AB, bh_BA ] }

	RECIPROCAL_BEST_HITS( recip_bh_ch )


	// Output for other workflows
	emit:
    rbh      = RECIPROCAL_BEST_HITS.out
}
