process CHECK_INPUTS {

    label 'local_job'

	input:
		path A_faa,  stageAs : "A/*"
		path B_faas, stageAs : "Bs/*"
        val type

	output:
        path "A/*.faa",  emit: "A_faa"
        path "Bs/*.faa", emit: "B_faas"

    script:
        """
        for fasta in $A_faa $B_faas
        do
            if [ ! -s \$fasta ]
            then
                echo "File \$fasta does not exist or is empty"
                exit 1
            fi

            # Translate if necessary
            if [[ $type == "NT" ]]
            then
                echo "Translating \$fasta"
                faTranslate.sh \$fasta
            fi
        done
        """
}
