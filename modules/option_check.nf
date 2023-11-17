process option_check {
	label 'default'
	scratch params.scratch
    errorStrategy 'terminate'
    input:
    output:
        val 'true', emit: readystate
shell:
'''
#!/bin/bash
ERROR=0
# FAM file
if [ ! -f "!{params.fam}" ]; then
    echo "Cannot open FAM file. Please fix the filename given to --fam. You also might want to check the documentation about mounting external paths:" >/dev/stderr
    ERROR=1
fi
# Covars
if [ "!{params.more_covars}" != "." ] && [ ! -f "!{params.more_covars}" ]; then
    echo "Cannot open covariate file.  Please fix the filename given to --fam. You also might want to check the documentation about mounting external paths:" >/dev/stderr
    ERROR=1
fi
if [ "!{params.more_covars}" != "." ] && [ -f "!{params.more_covars}" ] && [ "!{params.more_covars_cols}" == "" ]; then
    echo "A covariate file has been specified with --more_covars but no covar columns were given using --more_covars_cols." >/dev/stderr
    ERROR=1
fi
# test if selected covars are available in provided covars file
if [ "!{params.more_covars}" != "." ] && [ -f "!{params.more_covars}" ]; then
    
    head -n1 "!{params.more_covars}" | gawk 'BEGIN { split("'"!{params.more_covars_cols}"'",cols,",") } { for (c in cols) { found=0; for (i=3; i<=NF; i++) { if ($i == cols[c]) {found=1; break} }; if (found==0) exit 1 }}'
    if [[ $? == 1 ]]; then
        echo "At least one covar column specified in --more_covars_cols is not found in provided covariate file." >/dev/stderr
        ERROR=1
    fi
fi
# Trait type
case "!{params.trait}" in
    binary) ;;
    quantitative) ;;
    *) echo "Unsupported trait type. Please specify 'binary' or 'quantitative' for --trait. " >/dev/stderr
       ERROR=1
       ;;
esac
exit $ERROR
'''
}
