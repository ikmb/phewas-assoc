process saige_step1 {
	tag "${params.collection_name}_${phenoflag}"
	scratch params.scratch
	label 'saige'
	cache 'lenient'

	input:
		tuple path(bed), path(bim), path(fam), path(logfile), path(covars), path(covars_cols), val(meta), path(phenofile)
	output:
		tuple val(phenoflag), val(meta), path(bed), path(bim), path(fam), path(saigeoutput), path(saigevariance)
	when: meta.valid
	script:
		def step1input = params.regenie_step1_input ? "${params.regenie_step1_input}" : "tmp"
		def TRAIT_ARGS = ''
		TRAIT_ARGS =  	(params.trait == 'binary') ? '--traitType=binary' :
						(params.trait == 'quantitative') ? '--traitType=quantitative --invNormalize=TRUE' : ''
		saige_covars = phenofile.getSimpleName() + '_combined_covars_pheno.csv'
		saigeoutput = phenofile.getSimpleName() + '.rda'
		saigevariance = phenofile.getSimpleName() + '.varianceRatio.txt'
		phenoflag = phenofile.getSimpleName()
		def raretestoption = params.saige_test_rare ? "--isCateVarianceRatio=TRUE --cateVarRatioMaxMACVecInclude=10.5,20.5,30.5 --cateVarRatioMinMACVecExclude=1.5,10.5,20.5" : "" //NOT TESTED
	"""
	export R_LIBS_USER=/dev/null

	sed 's/^chr//' ${bim.baseName}.bim >tmp.bim

	ln -s ${bed} tmp.bed
	ln -s ${fam} tmp.fam

	Rscript ${projectDir}/bin/combine_col_pheno.R ${phenofile} ${covars} ${params.trait} ${saige_covars}

    Rscript /usr/local/bin/step1_fitNULLGLMM.R     \
        --plinkFile=$step1input  \
        --phenoFile=$saige_covars \
        --phenoCol=${phenofile.getSimpleName()} \
        --covarColList=\$(cat ${covars_cols}) \
        --sampleIDColinphenoFile=IID \
        $TRAIT_ARGS \
        --outputPrefix=./${phenofile.getSimpleName()} \
        --nThreads=${task.cpus}	\
        --IsOverwriteVarianceRatioFile=TRUE \
		$raretestoption
	"""
}

process saige_step2 {
	scratch params.scratch
	tag "${params.collection_name}_${phenoflag}"
	label 'saige'
	publishDir params.output, mode: 'copy'
	cache 'lenient'

	input:
		tuple val(phenoflag), val(meta), path(bed), path(bim), path(fam), path(saigeoutput), path(saigevariance)
	output:
		path(outprefix)
	when: meta.valid
	shell:
		outprefix = params.collection_name + '_' + phenoflag + '_saige.txt'// + params.test
		def raretestoption = params.saige_test_rare ? "---minMAF=0 --minMAC=0.5 --cateVarRatioMaxMACVecInclude=10.5,20.5,30.5 --cateVarRatioMinMACVecExclude=5.5,10.5,20.5" : "--minMAF=0 --minMAC=20"

		"""
		Rscript /usr/local/bin/step2_SPAtests.R \
			--bimFile=$bim \
			--bedFile=$bed \
			--famFile=$fam \
			--AlleleOrder=ref-first \
			--SAIGEOutputFile=${outprefix} \
			--GMMATmodelFile=${saigeoutput} \
			--varianceRatioFile=${saigevariance}	\
			--is_Firth_beta=TRUE    \
			--pCutoffforFirth=0.05 \
			--is_output_moreDetails=TRUE    \
			--LOCO=FALSE \
			$raretestoption
		"""
}