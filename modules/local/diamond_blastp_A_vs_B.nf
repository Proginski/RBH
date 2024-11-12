process DIAMOND_BLASTP_A_VS_B {

	// BLASTp a fasta A against a fasta B.

	label 'parallel_job'

	publishDir "${params.outdir}/blast_out", pattern: "*_BLASTp_*_*.out"

	input:
		path A_faa, stageAs: 'A/*'
		path B_faa, stageAs: 'B/*'

	output :
		tuple path("${A_faa}"), path("${B_faa}"), path("${file(A_faa).getName()}_BLASTp_${file(B_faa).getName()}_*.out")

	"""
	chmod -R +x ${projectDir}/bin

	# Setup
	echo "cpus = ${task.cpus}"
	timestamp=\$(date -d "today" +"%Y%m%d_%H%M")
	A_base=\$(basename ${A_faa})
	B_base=\$(basename ${B_faa})

	# Create the required databases :
	diamond makedb --in ${B_faa} -d ${B_faa} --threads ${task.cpus}


	# ### BLASTp ### A against B
	# blastp header = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
	diamond blastp  --query ${A_faa} --db ${B_faa} --threads ${task.cpus} \
	-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen \
	--evalue 0.001 --unal 1 -k 0 --max-hsps 0 \
	-o \${A_base}_BLASTp_\${B_base}_\${timestamp}.out


	# Add a last 'qcovs' column
	add_qcovs.sh \${A_base}_BLASTp_\${B_base}_\${timestamp}.out > \${A_base}_BLASTp_\${B_base}_\${timestamp}.out.tmp
	mv \${A_base}_BLASTp_\${B_base}_\${timestamp}.out.tmp \${A_base}_BLASTp_\${B_base}_\${timestamp}.out
	"""
}