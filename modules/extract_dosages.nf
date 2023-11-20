process extract_dosage {
    scratch params.scratch
    tag "${params.collection_name}.${chrom}"

    label 'extract_dosage'

    input:
    tuple file (gz), file(tbi), val(chrom), val(filetype)

    output:
    tuple file("${chrom}.PLINKdosage.map"), file("${chrom}.PLINKdosage.gz"), val(chrom)

shell:
'''
TARGET=!{chrom}.PLINKdosage

FIFO=$TARGET.fifo
rm -f $FIFO
mkfifo $FIFO

(<$FIFO tail -n +2 | gawk '{print $1, $3, 0, $2}' | uniq) >$TARGET.map &

# unpack vcf.gz only once
bgzip -c -d !{gz} \
    | gawk -f !{baseDir}/bin/extract_dosages.awk \
    | tee $FIFO \
    | gzip >$TARGET.gz

wait

rm -f $FIFO
'''
}