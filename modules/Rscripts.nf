process regenie_plot_manhattan {
	tag "${regenie_sumstats.getSimpleName()}"
	scratch params.scratch
	label 'rplotting'
	label 'short_run'
	publishDir params.output, mode: 'copy'
	input:
		path(regenie_sumstats)
	output:
		path('*')
	shell:
		'''
		Rscript !{baseDir}/bin/plot_manhattan_regenie.R !{regenie_sumstats} !{regenie_sumstats.getSimpleName()}
        '''
}

process regenie_plot_qq {
	tag "${regenie_sumstats.getSimpleName()}"
	scratch params.scratch
	label 'rplotting'
	label 'short_run'
	publishDir params.output, mode: 'copy'
	input:
		path(regenie_sumstats)
	output:
		path('*')
	shell:
		'''
		Rscript !{baseDir}/bin/plot_qq_regenie.R !{regenie_sumstats} !{regenie_sumstats.getSimpleName()}
        '''
}
