process plink_assoc {
    tag "${params.collection_name}.${chrom}"
    label 'plink'
    scratch params.scratch
    input:
    tuple file(dosagemap), file(dosage), val(chrom)
    each tuple file(pheno), file(cols)
    each file fam

    output:
    file "${chrom}.assoc.dosage"

shell:
'''
MEM=!{task.memory.toMega()-1000}

plink --fam !{fam} --map !{dosagemap} --dosage !{dosage} skip0=2 skip1=0 skip2=1 format=3 case-control-freqs \\
    --covar !{pheno} --covar-name $(cat !{cols}) \\
    --allow-no-sex --ci 0.95 \\
    --out !{chrom} --memory $MEM

'''
}