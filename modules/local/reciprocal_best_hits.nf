process RECIPROCAL_BEST_HITS{

    label 'local_job'

	// Builds a two cols tsv which is the join of two reciprocal blast best hits tsv files.

	publishDir "${params.outdir}/orthologs/", mode: 'copy'

    input:
    tuple val(A_raw), val(B_raw), path(bh_AB, stageAs : "bh_AB/*"), path(bh_BA, stageAs : "bh_BA/*")

    output:
    tuple val(A_raw), val(B_raw), path("${A_raw}_vs_${B_raw}_orthologs.tsv")

    """
    touch ${A_raw}_vs_${B_raw}_orthologs.tsv

	awk '
        BEGIN{FS=OFS="\t"}
        
        # First file
        FNR==NR { bh[\$1]=\$2 }
        
        # Second file
        FNR!=NR && bh[\$2]==\$1

	' $bh_BA $bh_AB > ${A_raw}_vs_${B_raw}_orthologs.tsv
    """
}