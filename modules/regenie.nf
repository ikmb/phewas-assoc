process phenofile_from_fam {
	tag "${params.collection_name}"
	scratch params.scratch
	label 'base'

	input:
		val(readystate)
		path(assocfam)
	output:
		path('phenotype.txt')

	shell:
	'''
	#gawk  'NR==1  {print "FID\tIID\tPhenotype"}{print "0\t"$1"_"$2"\t"$6}' !{assocfam} > phenotype.txt #-v 'OFS= '  -F '\t'
	#gawk  'NR==1  {print "FID\tIID\tPhenotype"}{print $1"\t"$2"\t"$6}' !{assocfam} > phenotype.txt #-v 'OFS= '  -F '\t'
	#gawk  'NR==1  {print "FID\tIID\tPhenotype"}{if($1!="0") {$2=$1"_"$2; $1="0";} print $1"\t"$2"\t"$6}' !{assocfam} > phenotype.txt #-v 'OFS= '  -F '\t'
	gawk  'NR==1  {print "FID\tIID\tphenotype"}{if($1!="0") {$2=$1"_"$2; $1="0";} if($6=="-9") {$6="NA";} print $1"\t"$2"\t"$6}' !{assocfam} > phenotype.txt #-v 'OFS= '  -F '\t'
	'''
}

process split_input_phenofile {
	tag "${params.collection_name}"
	scratch params.scratch
	label 'base'

	input:
		val(readystate)
		path(phenofile)
	output:
		path('*.txt')

	shell:
	'''
	#!/bin/bash

	# Set the input filename
	input_file="!{phenofile}"

	# Loop through each column in the input file, starting at column 3 (the first two columns are FID and IID)
	for ((i=3;i<=\$(head -1 \$input_file | awk '{print NF}');i++))
	do
		# Set the output filename based on the column name
		column_name=\$(head -1 $input_file | cut -d ' ' -f \$i)
		output_file="\${column_name}.txt"

		# Create the output file with headers
		#echo -e "FID\tIID\t${column_name}" > "$output_file"

		# Append the data to the output file
		cut -d ' ' -f 1,2,\$i "$input_file" > "$output_file"
	done
	rm !{phenofile}
	'''
}

process regenie_step1 {
	tag "${params.collection_name}_${phenofile.getSimpleName()}"
	scratch params.scratch
	label 'regenie'
	cache 'lenient'

	input:
		tuple path(bed), path(bim), path(fam), path(logfile), path(covars), path(covars_cols), val(meta), path(phenofile)
	output:
		tuple path('fit_bin_out_*.loco*'), path('fit_bin_out_pred.list'), path(covars), path(covars_cols), val(meta), path(phenofile)

	when: meta.valid
	script:
	def step1input = params.regenie_step1_input ? "${params.regenie_step1_input}" : "tmp"
	def catcovars = params.cat_covars ? "--catCovarList ${params.cat_covars}" : ""
	def TRAIT_ARGS = ''
		TRAIT_ARGS =  	(params.trait == 'binary') ? '--bt --cc12' :
						(params.trait == 'quantitative') ? '--qt --apply-rint' : ''
	def BUILT_ARGS = ''
		BUILT_ARGS =  	(params.build == '37') ? '--par-region hg37' :
						(params.build == '38') ? '--par-region hg38' : ''
	if (!params.regenie_step1_input) {  
		"""
		sed 's/^chr//' ${bim.baseName}.bim >tmp.bim

		ln -s ${bed} tmp.bed
		ln -s ${fam} tmp.fam

		export R_LIBS_USER=/dev/null

		regenie \
			--step 1 \
			--bed $step1input \
			--threads ${task.cpus} \
			--covarFile ${covars} \
			--covarCol \$(cat ${covars_cols}) \
			$catcovars \
			--phenoFile ${phenofile} \
			--use-relative-path \
			--bsize 100 \
			$TRAIT_ARGS \
			$BUILT_ARGS \
			--lowmem \
			--loocv	\
			--lowmem-prefix tmp_rg \
			${params.additional_regenie_parameter} \
			--out fit_bin_out \
			--gz
			"""
	} else {
		"""
		sed 's/^chr//' ${bim.baseName}.bim >tmp.bim
		ln -s ${bed} tmp.bed
		ln -s ${fam} tmp.fam


		sed 's/^chr//' ${params.regenie_step1_input}.bim >tmp_input.bim
		ln -s ${params.regenie_step1_input}.bed tmp_input.bed
		ln -s ${params.regenie_step1_input}.fam tmp_input.fam

		export R_LIBS_USER=/dev/null

		regenie \
			--step 1 \
			--bed tmp_input \
			--threads ${task.cpus} \
			--covarFile ${covars} \
			--covarCol \$(cat ${covars_cols}) \
			$catcovars \
			--phenoFile ${phenofile} \
			--use-relative-path \
			--bsize 100 \
			$TRAIT_ARGS \
			$BUILT_ARGS \
			--lowmem \
			--loocv	\
			--lowmem-prefix tmp_rg \
			${params.additional_regenie_parameter} \
			--out fit_bin_out \
			--gz
		"""
	}
}
//--phenoCol "Phenotype" \

process regenie_step2 {
	scratch params.scratch
	tag "${params.collection_name}_${phenofile.getSimpleName()}"
	label 'regenie'
	publishDir params.output, mode: 'copy'
	cache 'lenient'

	input:
		tuple path(bed), path(bim), path(fam), path(logfile), path(locofiles), path(predlist), path(covars), path(covars_cols), val(meta), path(phenofile)
	output:
		path("${params.collection_name}_regenie_${params.test}*"), emit: regenie_allout
		path("*.regenie.gz"), emit: sumstat
	when: meta.valid
	script:
		def catcovars = params.cat_covars ? "--catCovarList ${params.cat_covars}" : ""
		outprefix = params.collection_name + '_regenie_' + params.test

		def TRAIT_ARGS = ''
			TRAIT_ARGS =  	(params.trait == 'binary') ? '--bt --cc12' :
							(params.trait == 'quantitative') ? '--qt --apply-rint' : ''
		def BUILT_ARGS = ''
			BUILT_ARGS =  	(params.build == '37') ? '--par-region hg37' :
							(params.build == '19') ? '--par-region hg19' :
							(params.build == '38') ? '--par-region hg38' : ''

		def TEST_ARGS = ''
			TEST_ARGS = (params.test == 'firth') ? '--firth --approx' :
						(params.test == 'spa') ? '--spa' : ''

		"""
		sed 's/^chr//' ${bim.baseName}.bim >tmp.bim

		ln -s ${bed} tmp.bed
		ln -s ${fam} tmp.fam

		regenie \
			--step 2 \
			--bed tmp \
			--threads ${task.cpus} \
			--covarFile ${covars} \
			--covarCol \$(cat ${covars_cols}) \
			$catcovars \
			--phenoFile ${phenofile} \
			--bsize 200 \
			$TRAIT_ARGS \
			$TEST_ARGS \
			$BUILT_ARGS \
			--pThresh ${params.pthresh} \
			--loocv	\
			--pred ${predlist} \
			--out ${outprefix} \
			${params.additional_regenie_parameter} \
			--gz
		"""
}
//#--covarCol PC{1:!{params.pca_dims}} \
process awk_regenie {
	tag "${regenie_sumstats.getSimpleName()}"
	scratch params.scratch
	label 'base'
	publishDir params.output, mode: 'copy'
	input:
		path(regenie_sumstats)
	output:
		path(regenie_newname)
	shell:
		regenie_newname = regenie_sumstats.getSimpleName() + '_plainpvalue.regenie.gz'
		'''
		zcat !{regenie_sumstats} | gawk 'NR==1{$14="p.value"; $4="ALLELE1"; $5="ALLELE2"; $6="A2FREQ"; print $0}NR>1{$14=10**-$12; print $0}' |  gawk 'NR>1 {if($14>1) next} 1' | gzip > !{regenie_newname}
		'''
}