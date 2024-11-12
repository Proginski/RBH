process FASOMERECORDS_B_BH{

    // Turns the second column of a *best_hits.tsv into a list to subset fasta B (target)

    input : 
    tuple path(A_faa, stageAs : "A/*"), path(B_faa, stageAs : "B/*"), path(bh)

    output :
    tuple path(A_faa), path(B_faa), path("${B_faa}.bh.faa")

    """
    awk -F"\t" '{print \$2}' $bh > ${B_faa}.bh.txt

    faSomeRecords $B_faa ${B_faa}.bh.txt ${B_faa}.bh.faa
    """
}