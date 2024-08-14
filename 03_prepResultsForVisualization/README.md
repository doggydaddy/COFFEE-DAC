# Converting permutation results to edges to plot for the visualizer

First step is the take the results from the permutation tests and convert them into a format for the visualizer.

## Create index dictionary

Firstly, we must create dictionary to convert permutation test result indices to ijk/xyz indices. 

Given the brain mask you used to during preprocessing/analysis (*\<brain_mask.nii\>*), create conversion look-up file. This is easiest done using the AFNI program:

        3dmaskdump -mask <brain_mask.nii> -index -xyz -o <index_dict.txt> <brain_mask.nii>

## Formatting permutation results

The visualizer takes a \*.csv file with edges to visualize (filtered). There are several ways to accomplish this.

### AFNI routines

If you used AFNI programs (e.g. *3dttest++*), the results will be presented in a
single (enormous) \*.nii file (3D+t, where *t* is the number of voxels). This
would be incredibly large to to store and is typically clipped so that each *t*
is stored in a separate file.

Several help routines are available to convert this result format into a \*.csv
file for the visualizer:

        ./find_non_empty_datasets.sh <data directory> <optional: xyz/ijk>

locates all volumes in a directory with non-zero values, and using the *AFNI*
program *3dmaskdump* to convert these volumes into text files with coordinates
(optional parameter, xyz/ijk, with xyz being the default if not specified), puts
the filenames in a file list moves the converted volumes into a intermediate
directory.

A *R* script then can be called to convert the non-zero volume conversions in
the intermediate directory to a single \*.csv file ready for loading.

        ./construct_edge_table.Rscript

These scripts require *R* and and *AFNI* installed, which should not be an issue
since these scripts are only used when statistical results are in *AFNI* format.

### cudaPerm results

