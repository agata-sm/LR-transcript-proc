 #!/usr/bin/env bash

let no_args=$#


declare -r -A array1


for var in "$@"
do

    mod=${var//[\[,\]]/}

    echo "$mod"

    array1+=("$mod")
done



declare -r -A array2

idx=0

while [ "$idx" -lt "$no_args" ]
do

    path=${array1[$idx]}
    smpl=${array1[$(($idx+1))]}

	printf '%s\t%s\n' "${path}" "${smpl}" 

    let idx=$idx+2

done > "sample_sheet.tsv"


