
process BLAST_BEST_HITS {

	// Returns a two cols tsv files (query, best_hit) from a blast output, with a filter on qcovs (col 14).
	
	publishDir "${params.outdir}/best_hits", pattern: "*.best_hits.tsv"

	input:
		tuple path(A_faa, stageAs : "A/*"), path(B_faa, stageAs : "B/*"), path(blast_out)
		
	output:
		tuple path(A_faa), path(B_faa), path("${blast_out}.best_hits.tsv")
		
	"""
	# The following is based on the principle that blast outputs are sorted for each query by decreasing e-values
	awk '
		BEGIN {FS=OFS="\t"}
		
		# For HSP (non-commment) lines,
		# with a query coverage by subject >= 70%, and an unprinted query,
		# Save the query and print the query-subject pair.
		\$0 !~ /^#/ && \$14 >= 70 && (!(\$1 in data )) { data[\$1]=1 ; print \$1,\$2 }
		
	' $blast_out > ${blast_out}.best_hits.tsv
	"""
}