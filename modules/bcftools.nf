process gen_r2_list {
	label 'bcftools'
    label 'short_run'
	scratch params.scratch
    tag "${params.collection_name}.${chrom}"

    input:
        tuple path(vcf), path(tbi), val(chrom)

    output:
        path("r2-include.${chrom}")

	script:
    def bcftoolfilter = params.additional_bcftools_arg ? "${params.additional_bcftools_arg}" : ""

        """
        set +e
        bcftools view $bcftoolfilter ${vcf} | bcftools query -i '${params.null_filter}' -f '%CHROM:%POS:%REF:%ALT\\n' >r2-include.${chrom}
        if [ \$? -ne 0 ]; then
            echo Filter not found or genotyped-only data available.
            bcftools view $bcftoolfilter ${vcf} | bcftools query -f '%CHROM:%POS:%REF:%ALT\\n' >r2-include.${chrom}
        fi
        exit 0
        """
}

process split_vcf_ranges {
    tag "${params.collection_name}.${chrom}"
	label 'bcftools'
    label 'short_run'
	scratch params.scratch
    input:
        tuple file(vcf), file(tbi), val(chrom)
    output:
        tuple file(vcf), file(tbi), file(field), val(chrom), file("0*")

    shell:
        '''
        bcftools query -f "%POS\\n" !{vcf} >positions
        bcftools view -H !{vcf} > viewed
        head viewed -n1 | cut -f9 >field
        split -d -a 8 -l !{params.saige_chunk_size} positions '0' #-d --suffix-length=8 --lines=!{params.saige_chunk_size}
        '''
}

process bcftoolsfilter {
    tag "${params.collection_name}.${chrom}"
	label 'bcftools'
    label 'short_run'
	scratch params.scratch
    input:
        tuple file(vcf), file(tbi), val(chrom)
    output:
        tuple file(outputname), file(outputname_tbi), val(chrom)

    script:
        def bcftoolfilter = params.additional_bcftools_arg ? "${params.additional_bcftools_arg}" : ""
        outputname = chrom +  '_bcffiltered.vcf.gz'
        outputname_tbi = chrom +  '_bcffiltered.vcf.gz.tbi'
        """
        bcftools view $bcftoolfilter ${vcf} -Oz --threads ${task.cpus} -o $outputname
        bcftools index $outputname --tbi --threads ${task.cpus}
        """
}