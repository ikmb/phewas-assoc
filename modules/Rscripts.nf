/*
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
*/
process plot_manhattan_qq {
	tag "${assoc_tool}_${sumstats.getSimpleName()}"
	scratch params.scratch
	label 'rplotting'
	label 'short_run'
	publishDir params.output, mode: 'copy'
	input:
		tuple val(assoc_tool), path(sumstats)
	output:
		path('*')
	script:
		def output_name = params.manhattan_output_name ? "${sumstats.getSimpleName().replaceAll('_plainpvalue','')}_${params.manhattan_output_name}" : "${sumstats.getSimpleName().replaceAll('_plainpvalue','')}"
		def title_text = params.manhattan_title_text ? params.manhattan_title_text : "${assoc_tool.toUpperCase()} association results"
		"""
		Rscript ${baseDir}/bin/make_manhattan_plot.R ${assoc_tool} \
													${sumstats} \
													"${title_text}" \
													${output_name}
        """
}