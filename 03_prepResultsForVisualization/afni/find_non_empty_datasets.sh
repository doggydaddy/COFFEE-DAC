#!/bin/bash

data_dir=$1
mask_file=$2

if [[ $# -eq 0 ]]; then
    echo 'USAGE ./find_non_empty_dataset.sh <data directory> <mask file> <ijk/xyz>'
fi

if [ $# -gt 2 ]; then
    ijk_or_xyz=$3
else
    echo 'assuming xyz'
    ijk_or_xyz='xyz'
fi

# in case data folder does not only contain ttest .niis
all_files=`find $data_dir -name 'xROI*.nii'`

# save a list of .niis with non-zero values
touch non_zero_vols.txt

for f in $all_files
do
    admax=`3dinfo -dmax $f`
    if (( $admax > 0 )); then
        echo $f >> non_zero_vols.txt  
        if [ "$ijk_or_xyz" = "ijk" ]; then
            3dmaskdump -mask $mask_file -o ${f%.nii}.dump $f
        elif [ "$ijk_or_xyz" = "xyz" ]; then
            3dmaskdump -mask $mask_file -xyz -o ${f%.nii}.dump $f
        else
            echo 'either ijk or xyz required as third argument'
        fi
    fi
done

mv data/*.dump intermediates/

