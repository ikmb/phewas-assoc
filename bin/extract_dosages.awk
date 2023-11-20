BEGIN{
    printf "CHR BP SNP A1 A2 INFO"
    gpfield = 0
    gtfield = 0
}
{
    # skip header lines
    if (($1 ~ /^##/)) {
    } else {
        # Copy header line
        if (($1 ~ /^#CHROM/)) {
            # Process one FID_IID-part into "0 FID_IID", so Plink thinks there
            # are really two values
            for (i=10; i<=NF; i++)
                printf " "0" "$i
            printf "\\n"
        } else {
            # Find field index of GP and GT
            n = split($9, fields, ":")
            for(i=1; i in fields; i++) {
                if(fields[i]=="GP") {
                    gpfield = i
                }
                if(fields[i]=="GT") {
                    gtfield = i
                }
            }

            if(gtfield == 0 && gpfield == 0) {
                print "Error: Cannot find a GT or GP field in input VCF." >"/dev/stderr"
                exit 1
            }

            # Extract r^2 score from INFO field, set INFO to 1.0 if none found
            {
                n = split($8, fields, /(;|=)/)
                info = "R2=1.0" # default to 1 if no usable info score is found
                for(i=1; i in fields; i++) {
                    if(fields[i] == "INFO" || fields[i] == "DR2" || fields[i] == "R2") {
                        info = fields[i+1]
                        break
                    }
                }
            }

            # keep a dictionary of variants to skip duplicates
            if(!seenbefore[$1":"$2":"$4":"$5]) {
                printf $1" "$2" "$1":"$2":"$4":"$5" "$4" "$5" "info

                # split value field
                for (i=10; i<=NF; i++) {
                    n=split($i,array,":")
                    # Choose GP over GT
                    if(gpfield == 0) {
                        n=split($i,array,/[\\/|\\|]/)
                        if(n==2) {
                            # "convert" hard calls to probabilities
                            if(array[1]+array[2]==0) { printf " 1 0 0" }
                            else if(array[1]+array[2]==1) { printf " 0 1 0" }
                            else { printf " 0 0 1" }
                        
                        } else if(n==1) {
                            # male samples on chrX, only one genotype
                            printf " "$i " 0 "$i
                        } else {
                            print "VCF Error: GP or GT field not present or malformatted" >"/dev/stderr"
                            exit 1
                        }
                    } else {

                        # Missing fields are specified as ".", not splittable by ","
                        if(array[gpfield] == ".") {
                            printf " . . ."
                        } else {
                            n=split(array[gpfield],array_GP,",")
                            if(n==3) {
                                printf " "array_GP[1]" "array_GP[2]" "array_GP[3]
                            } else {
                                # males on chrX
                                printf " "array_GP[1]" 0 "array_GP[2]
                            }

                        }
                    }
                }

                printf "\\n"

                seenbefore[$1":"$2":"$4":"$5] = 1
            }
        }

    }
}