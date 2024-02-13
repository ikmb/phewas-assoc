#!/bin/bash

# Check if the input file is provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file="$1"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "File '$input_file' not found."
    exit 1
fi

# Check each column except the first and second
gawk 'BEGIN { FS = " "; invalid = 0 }
     NR > 1 {
        for (i = 3; i <= NF; i++) {
            if ($i != "1" && $i != "2" && $i != "NA") {
                invalid = 1
                break
            }
        }
     }
     END {
        if (invalid == 1) {
            print "Invalid values found in the input file. If trait is set to binary valid values are 1/2/NA"
            exit 1
        } else {
            print "All values in the input file are valid."
            exit 0
        }
     }' "$input_file"