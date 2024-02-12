process option_check {
	label 'base'
	scratch params.scratch
    errorStrategy 'terminate'
    input:
    output:
        val 'true', emit: readystate
shell:
'''
#!/bin/bash
# test if selected covars are available in provided covars file
if [ "!{params.more_covars}" != "." ] && [ -f "!{params.more_covars}" ]; then
    head -n1 "!{params.more_covars}" | gawk 'BEGIN { split("'"!{params.more_covars_cols}"'",cols,",") } { for (c in cols) { found=0; for (i=3; i<=NF; i++) { if ($i == cols[c]) {found=1; break} }; if (found==0) exit 1 }}'
    if [[ $? == 1 ]]; then
        echo "At least one covar column specified in --more_covars_cols is not found in provided covariate file." >/dev/stderr
        exit 1
    fi
fi
# Trait type
case "!{params.trait}" in
    binary) ;;
    quantitative) ;;
    *) echo "Unsupported trait type. Please specify 'binary' or 'quantitative' for --trait. " >/dev/stderr
        exit 1
        ;;
esac
# test if phenofile is only containing 1/2/NA values and is space separated should it be used for binary traits.
if [ "!{params.trait}" == "binary" ] && [ "!{params.phenofile}" != "false" ]; then
    if gawk '$1=$1' OFS='\\t' "!{params.phenofile}" | grep -qP '\\S+\\s+\\S+'; then
        # echo "File "!{params.phenofile}" is space-separated."

        bash !{projectDir}/bin/test_binary_phenofile.sh "!{params.phenofile}"

    else
        echo "File "!{params.phenofile}" is not space-separated, exiting now!" >/dev/stderr
        exit 1
    fi
fi
'''
}
