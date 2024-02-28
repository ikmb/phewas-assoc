process plink2_assoc {
    tag "${params.collection_name}_${phenotype}_${chrom}"
	scratch params.scratch
    label 'plink2'

    input:
    tuple val(meta), path(phenofile), path(vcf), path(tbi), val(chrom), val(filetype), path(covars), path(covars_cols)
    each path(fam)

    output:
        path('*'), emit: all
        tuple val(phenotype), path('*.glm*'), emit: plinksumstats
    when: 
        meta.valid
    script:
        phenotype = phenofile.getSimpleName()
        output_name = chrom + '.plink2_assoc_' + phenotype
        //allow-no-covars 
        def glmoptions = params.plink2_glm_options ? "--glm ${params.plink2_glm_options}" : "--glm omit-ref hide-covar --mac 20"
        def memory = task.memory.toMega()-1000
    """
        #Create a fam file to update sex information
        echo "#FID IID PAT MAT SEX PHE" >new-fam
        awk '{\$3="0";\$4="0";if(\$1!="0") {\$2=\$1"_"\$2; \$1="0";} print \$0}' ${fam} >>new-fam

        plink2 --vcf ${vcf} \
            --threads ${task.cpus} \
            --memory $memory \
            --update-sex new-fam \
            --out ${output_name} \
            --chr ${chrom} \
            --allow-extra-chr \
            $glmoptions \
            --pheno-name ${phenotype} \
            --covar ${params.collection_name}.covars \
            --covar-name \$(cat ${covars_cols}) \
            --split-par "b${params.build}" \
            --pheno ${phenofile}
    """
}
//TODO: include dosage like so: --vcf 6.ap_prf.vcf.gz dosage=GP

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
        merged_sumstats = phenotype + '_plink2_glm_sumstats.tsv.gz'
        
        '''
        file_list=!{sumstats.join(',')}

        header=$(head -n 1 !{sumstats[0]} | sed '/^\s*$/d')  

        # Create the output file with the header
        echo "$header" > tmp.tsv

        #Cat in all files without their header line
        for i in ${file_list//,/ }
        do
            awk 'NR > 1' "$i" >> tmp2.tsv
        done

        #sort that all bp are in correct order, first chromosome sorting
        awk 'NR<2{print $0;next}{print $0| "sort -k1,1n -k2,2n"}' tmp2.tsv 

        #now split results into two files should there be chromosomes like "X" or "Y"
        sed -e '/^[\\t\\v\\f ]*[^0-9]/ d' tmp2.tsv  | sort -k1,1n -k2,2n > numbers.tsv
        sed -e '/^[\\t\\v\\f ]*[0-9]/ d' tmp2.tsv  | sort -k1,1d -k2,2n > others.tsv

        #now merge first header, then numerical chromosomes and then alphabetical chromosomes
        cat tmp.tsv numbers.tsv others.tsv | gzip > !{merged_sumstats}

        '''
}
