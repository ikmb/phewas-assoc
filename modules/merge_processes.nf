process merge_plink_results {
    tag "${params.collection_name}"
	label 'base'
	scratch params.scratch
    publishDir params.output, mode: 'copy'

    input:
    path(stats)// from for_plink_results.collect()

    output:
    path("${params.collection_name}.Plink.stats") //into for_lift_sumstats_plink

shell:
'''
# extract first line, convert tabs to space
head -n1 !{stats[0]} | tr -s '\t ' ' ' | xargs >!{params.collection_name}.Plink.stats
ls !{stats} | sort -n | xargs -n1 tail -n +2 | gawk '{if(substr($1,1,3)!="chr"){$1="chr"$1} $2=$1":"$3":"$4":"$5; print}'>>!{params.collection_name}.Plink.stats
'''
}

process merge_r2 {
    tag "${params.collection_name}"
	scratch params.scratch
	label 'base'

    input:
    path(r2)

    output:
    path("r2-include.sorted")

	shell:
	'''
		cat r2-include.* >r2-include
		export TMPDIR=.
		gawk '$0 !~ /^chr/ {$1="chr"$1} {print}' r2-include | sort >r2-include.sorted
	'''
}
