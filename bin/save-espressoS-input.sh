 #!/usr/bin/env bash

echo "cmdline args"
echo "${[@]}"


let no_args=$#


declare -a array1


for var in "$@"
do

    mod=${var//[\[,\]]/}

    echo "$mod"

    array1+=("$mod")
done



idx=0

while [ "$idx" -lt "$no_args" ]
do

    path=${array1[$idx]}
    smpl=${array1[$(($idx+1))]}

	printf '%s\t%s\n' "${path}" "${smpl}" 

    let idx=$idx+2

done > "sample_sheet.tsv"


