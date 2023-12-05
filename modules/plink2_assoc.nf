process plink2_assoc {
    tag "${params.collection_name}_${phenotype}"
	scratch params.scratch
    //publishDir params.output, mode: 'copy'
    label 'plink2'

    input:
    tuple val(meta), path(phenofile), path(vcf), path(tbi), val(chrom), val(filetype), path(covars), path(covars_cols)

    output:
        path('*'), emit: all
        tuple val(phenotype), path('*.glm*'), emit: plinksumstats
    when: meta.valid
    script:
        phenotype = phenofile.getSimpleName()
        output_name = chrom + '.plink2_assoc_' + phenotype
        def glmoptions = params.plink2_glm_options ? "--glm ${params.plink2_glm_options}" : "--glm omit-ref hide-covar"
        def memory = task.memory.toMega()-1000
    """
        plink2 --vcf ${vcf} dosage=GP \
            --threads ${task.cpus} \
            --memory $memory \
            --out ${output_name} \
            $glmoptions \
            --pheno-name ${phenotype} \
            --covar ${params.collection_name}.covars \
            --covar-name \$(cat ${covars_cols}) \
            --pheno ${phenofile}
    """
}
//MEM=${task.memory.toMega()-1000}
process plink2_assoc_merge {
    tag "${params.collection_name}_${phenotype}"
	scratch params.scratch
    label 'base'
    publishDir params.output, mode: 'copy'
    
    input:
        tuple val(phenotype), path(sumstats)
        
    output:
        tuple val(phenotype), path(merged_sumstats)
    shell:
        merged_sumstats = phenotype + '_plink2_glm_sumstats.tbl'
        
        '''
        gawk 'BEGIN { firstfile = 1 }
        FNR == 1 && firstfile == 1 {
            header = $0
            firstfile = 0
            print header
        }
        FNR > 1 || FILENAME ~ /\\.glm$/ {
            print $0
        }' *.glm* | awk 'NR<2{print $0;next}{print $0| "sort -k1,1n -k2,2n"}' > !{merged_sumstats}
        '''
}
