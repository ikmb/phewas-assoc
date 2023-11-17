process make_covars {
	scratch params.scratch
    publishDir params.output, mode: 'copy'
	label 'base'
    input:

        path(evec)
        path(inc_fam)

    output:
        tuple path("${params.collection_name}.covars"), path("${params.collection_name}.covar_cols"), emit: covars

    shell:
    '''
    # Re-format sample ID and family ID to VCF rules
    gawk '{if($1!="0") {$2=$1"_"$2; $1="0";} print $0}' !{inc_fam}  | sort >new-fam

    gawk 'FNR==NR{samples[$2];next} {if($2 in samples) { for(i=1;i<=(12);i++) {printf "%s%s", $i, (i<12?OFS:ORS)}}}' new-fam !{evec} >filtered-evec
    if [ -f "!{params.more_covars}" ]; then
        gawk 'FNR==NR{samples[$2];next} {if($2 in samples) { for(i=1;i<=(NF);i++) {printf "%s%s", $i, (i<NF?OFS:ORS)}}}' new-fam !{params.more_covars} | sort >filtered-covars
        
        # identify column indices of selected covars
        # NOTE: We checked above that all selected columns are available, so the list won't be empty
        COLSTR=$(head -n1 "!{params.more_covars}"  | gawk 'BEGIN { split("'"!{params.more_covars_cols}"'",cols,","); colstr="" } { for (c in cols) { for (i=3; i<=NF; i++) { if ($i == cols[c]) {colstr=colstr "," i} }}} END { print substr(colstr,2) }')
        
        # extract header
        cut -f$COLSTR -d" " !{params.more_covars} | awk 'NR<=1' >covars-column
        # fill values
        cut -f$COLSTR -d" " filtered-covars >>covars-column
    fi
    EVEC_LINES=$(wc -l <filtered-evec)
    FAM_LINES=$(wc -l <new-fam)
    if [ $EVEC_LINES -ne $FAM_LINES ]; then
        echo I could not find covariates in !{evec} for every sample specified in !{inc_fam}. Please check.
        echo Aborting.
        exit 1
    fi
    echo "FID IID" PC{1..!{params.pca_dims}} >evec.double-id.withheader

    cat filtered-evec >>evec.double-id.withheader


    # Merge both, replace space runs with single tabs for SAIGE
    touch covars-column
    if [[ -s "!{params.more_covars}" ]]; then
        echo PC{1..!{params.pca_dims}}, !{params.more_covars_cols} | sed 's/\\ ,/\\,/g' | tr ' ' ,  | sed 's/\\,\\,/\\,/g' >!{params.collection_name}.covar_cols
    else
        echo PC{1..!{params.pca_dims}} | tr ' ' , >!{params.collection_name}.covar_cols
    fi

    header=$(head -n 1 covars-column)
    if [[ "$header" == "FID\tIID"* ]]; then
        awk '{$1="";$2="";print $0}' covars-column | sed 's/  //g' > covars-column_temp
        paste -d" " evec.double-id.withheader covars-column_temp | tr -s ' ' \\\\t > !{params.collection_name}.covars > !{params.collection_name}.covars
    else
        paste -d" " evec.double-id.withheader covars-column | tr -s ' ' \\\\t >!{params.collection_name}.covars
    fi
    '''
}

process make_covars_nopca {
	scratch params.scratch
    publishDir params.output, mode: 'copy'
	label 'base'
    input:
        val(readystate)
        path(fam)

    output:
        tuple path("${params.collection_name}.covars"), path("${params.collection_name}.covar_cols"), emit: covars

    shell:
    '''
    # Re-format sample ID and family ID to VCF rules
    gawk '{if($1!="0") {$2=$1"_"$2; $1="0";} print $0}' !{fam}  | sort >new-fam

    if [ -f "!{params.more_covars}" ]; then
        gawk 'FNR==NR{samples[$2];next} {if($2 in samples) { for(i=1;i<=(NF);i++) {printf "%s%s", $i, (i<NF?OFS:ORS)}}}' new-fam !{params.more_covars} | sort >filtered-covars
        
        # identify column indices of selected covars
        # NOTE: We checked above that all selected columns are available, so the list won't be empty
        COLSTR=$(head -n1 "!{params.more_covars}"  | gawk 'BEGIN { split("'"!{params.more_covars_cols}"'",cols,","); colstr="" } { for (c in cols) { for (i=3; i<=NF; i++) { if ($i == cols[c]) {colstr=colstr "," i} }}} END { print substr(colstr,2) }')
        
        # extract header
        cut -f$COLSTR -d" " !{params.more_covars} | awk 'NR<=1' >covars-column
        # fill values
        cut -f$COLSTR -d" " filtered-covars >>covars-column
    fi
    touch covars-column

    echo "FID IID" >evec.double-id.withheader

    gawk '{print $1" "$2}' new-fam >>evec.double-id.withheader


    # Merge both, replace space runs with single tabs for regenie
    
    if [[ -s "!{params.more_covars}" ]]; then
        echo !{params.more_covars_cols} | sed 's/\\ ,/\\,/g' | tr ' ' ,  | sed 's/\\,\\,/\\,/g' >!{params.collection_name}.covar_cols
    else
        echo '' | tr ' ' , >!{params.collection_name}.covar_cols
    fi

    header=$(head -n 1 covars-column)
    if [[ "$header" == "FID\tIID"* ]]; then
        gawk '{$1="";$2="";print $0}' covars-column | sed 's/  //g' > covars-column_temp
        paste -d" " evec.double-id.withheader covars-column_temp | tr -s ' ' \\\\t > !{params.collection_name}.covars > !{params.collection_name}.covars
    else
        paste -d" " evec.double-id.withheader covars-column | tr -s ' ' \\\\t >!{params.collection_name}.covars
    fi
    '''
}
