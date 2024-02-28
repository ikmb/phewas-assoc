process make_covars {
	scratch params.scratch
    publishDir params.output, mode: 'copy'
	label 'tidyverse'
    input:

        path(evec)
        path(inc_fam)

    output:
        tuple path("${params.collection_name}.covars"), path("${params.collection_name}.covar_cols"), emit: covars

    script:
        def withcovars = (params.more_covars && params.more_covars != ".") ? true : false
    """
        Rscript ${baseDir}/bin/make_covars.R $inc_fam ${params.more_covars} "${params.more_covars_cols}" $evec ${params.pca_dims} $withcovars ${params.collection_name}
    """
}

process make_covars_nopca {
	scratch params.scratch
    publishDir params.output, mode: 'copy'
	label 'tidyverse'
    input:
        val(readystate)
        path(inc_fam)

    output:
        tuple path("${params.collection_name}.covars"), path("${params.collection_name}.covar_cols"), emit: covars

script:
        def withcovars = (params.more_covars && params.more_covars != ".") ? true : false
    """
        Rscript ${baseDir}/bin/make_covars.R $inc_fam ${params.more_covars} "${params.more_covars_cols}" "nopca" ${params.pca_dims} $withcovars ${params.collection_name}
    """
}