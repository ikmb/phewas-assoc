process prune {

    label 'plink2'
	scratch params.scratch
    publishDir params.output, mode: 'copy'

    input:
        tuple path(bed), path(bim), path(fam), path(logfile)
        path(r2)
		path("include-r2-variants")

    output:
        tuple file("${params.collection_name}.pruned.bed"), file("${params.collection_name}.pruned.bim"), file("${params.collection_name}.pruned.fam"), file("${params.collection_name}.pruned.log")// into for_saige_null, for_liftover_pruned, for_generate_pcs

	shell:
    '''
    echo Generating PCA SNP List file for variant selection
    R2=!{r2}
    MEM=!{task.memory.toMega()-1000}
    if [ ! -s $R2 ]; then
        # No entries in include list. This means that either we have genotyped-only
        # data where no R2 score is present or the filter is malformed. We're assuming
        # genotyped-only data here. In this case, populate the r2 list.
        cut -f2 -d" " !{bim} | sort >new-r2
        R2=new-r2
    fi
    plink2 --memory $MEM --threads !{task.cpus} --bed !{bed} --bim !{bim} --fam !{fam} --extract $R2 --indep-pairwise !{params.pruning_parameter} --out _prune
    plink2 --memory $MEM --threads !{task.cpus} --bed !{bed} --bim !{bim} --fam !{fam} --extract _prune.prune.in --maf 0.05 --make-bed --out intermediate
    plink2 --memory $MEM --threads !{task.cpus} --bfile intermediate --chr 1-22 --extract include-r2-variants --output-chr chrM --make-bed --out !{params.collection_name}.pruned
    rm intermediate*
'''
}

process make_plink {
    label "plink2"

    scratch params.scratch
    tag "${params.collection_name}.$chrom"

    input:
        tuple path(vcf), path(tbi), val(chrom), val(filetype)
        each path(fam)

    output:
        tuple path("${params.collection_name}.${chrom}.bed"), path("${params.collection_name}.${chrom}.bim"), path("${params.collection_name}.${chrom}.fam"), path("${params.collection_name}.${chrom}.log"),val(chrom) //into for_liftover, for_merge_b38

    script:
        def MEM=task.memory.toMega()-1000
        if (params.fam) {  
            """
            # Generate double-id FAM
            

            #this was gawk originally:
            awk '{\$3="0";\$4="0";if(\$1!="0") {\$2=\$1"_"\$2; \$1="0";} print \$0}' ${fam} >new-fam
            #awk '{if(\$1!="0") {\$1="0";} print \$0}' ${fam} >new-fam

            plink2  --vcf ${vcf} \
                    --const-fid 0 \
                    --memory $MEM \
                    --max-alleles 2 \
                    --keep-nosex \
                    --chr ${chrom} \
                    --allow-extra-chr \
                    --pheno new-fam \
                    --pheno-col-nums 4 \
                    --update-sex new-fam col-num=3 \
                    --split-par "b${params.build}" \
                    --output-chr chrM \
                    --make-bed \
                    --out ${params.collection_name}.${chrom}

            mv ${params.collection_name}.${chrom}.bim old_bim

            #this was gawk originally:
            awk '\$1 !~ /^chr/ {\$1="chr"\$1} {\$2=\$1":"\$4":"\$6":"\$5; print}' <old_bim >${params.collection_name}.${chrom}.bim
            # Might need some "chr" prefixing here
            """
        }else{
            """
            # Generate double-id FAM
            #this was gawk originally:
            awk '{if(\$1!="0") {\$2=\$1"_"\$2; \$1="0";} print \$1" "\$2" "\$3}' ${fam} >new-fam
            #awk '{if(\$1!="0") {\$1="0";} print \$0}' !{fam} >new-fam

            plink2  --vcf ${vcf} \
                    --const-fid 0 \
                    --memory $MEM \
                    --max-alleles 2 \
                    --keep-nosex \
                    --chr ${chrom} \
                    --allow-extra-chr \
                    --pheno new-fam \
                    --pheno-col-nums 1 \
                    --split-par "b${params.build}" \
                    --output-chr chrM \
                    --make-bed \
                    --out ${params.collection_name}.${chrom}

            mv ${params.collection_name}.${chrom}.bim old_bim

            #this was gawk originally:
            awk '\$1 !~ /^chr/ {\$1="chr"\$1} {\$2=\$1":"\$4":"\$6":"\$5; print}' <old_bim >${params.collection_name}.${chrom}.bim
            # Might need some "chr" prefixing here
            """
        }
}

process merge_plink {
    tag "${params.collection_name}"
	scratch params.scratch
    publishDir params.output, mode: 'copy'
    label 'plink2'

    input:
    file(filelist)

    output:
    tuple file("${params.collection_name}.bed"), file("${params.collection_name}.bim"), file("${params.collection_name}.fam"), file("${params.collection_name}.log")// into for_prune //, for_liftover


    shell:
	'''
		ls *.bed | xargs -i -- basename {} .bed | tail -n +2 >merge-list
		MEM=!{task.memory.toMega()-1000}
		FIRSTNAME=$(ls *.bed | xargs -i -- basename {} .bed | tail -n +1 | head -n 1)
        
        plink2  --bfile $FIRSTNAME \
                --threads !{task.cpus} \
                --memory $MEM \
                --pmerge-list merge-list bfile \
                --make-bed \
                --keep-nosex \
                --indiv-sort none \
                --output-chr 26 \
                --chr 1-22, X, MT, par1, par2 \
                --out !{params.collection_name}
        #trying to remake the plink1.9 output with parameter --output-chr 26
        #remove CHR Y as regenie cannot handle it
    '''
}
// TODO: ADD "--sort-vars", needs a pgen output first, that needs then to be sorted and only then output with --make-bed