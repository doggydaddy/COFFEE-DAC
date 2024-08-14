# CUDA accelerated permutation testing

This subroutines performs connection-wise permutation tests.

The statistical tests are performed independently for each voxel, so the tests
can be greatly accelerated using GPU.

# Details

For each test, all N subjects included in the tests must be loaded into memory.
For easy memory management, prior to performing the tests, the permutations are
generated *a-priori*.

## Generating permutations

The permutations are generated using a python program *createPerm.py*. It
generates one-hot labels of $n$ permutations of nA and nB subjects, where $nA$
and $nB$ is the number of subjects in group A and group B respectively.

Calling the program is simple:

        python createPerm.py <nr. permutations> <nA> <nB> <output txt file>

**Note** that the original statistical test groupings is NOT known to the
program. So the first row of the generated text file must be replaced with the
original labels, otherwise the test will not be as intended! The ordering is
from left-to-right from the top-to-down order of appearance in *filelist.txt*

## Performing permutation tests

### Parsing file list

The file list is a simple list of files to be processed in a single column with
one file name per row. The order does matter as the intended statistical test
groupings must be in order of the *first* row in the permutations file. 

For example, the filelist may be initiated using:

        find . -name 'groupA_subj*_connectivityMatrix.txt' | sort > filelist.txt
        find . -name 'groupB_subj*_connectivityMatrix.txt' | sort >> filelist.txt

The program *parseFileList* parses all subjects in the *filelist* from *start
index* to *end index*, and outputs it into a separate *output file*.

        parseFileList <filelist> <start index> <end index> <output file>

Each row contains the subject voxel data from start to end index, in order
left-to-right according to the input file lists top-to-down order.

# Simulated data for testing

The iPython notebook *genTestData.ipynb* can be used to generate simulated data
for testing purposes. Simulated data with their names can be found and changed
in *subject_list*. *N* is the number of simulated voxels (not connections).
*subject_list_A* contains a sub-list, where each subject have certain
connections artificially truncated to above a threshold (0.8). Otherwise all
simulated connection values are uniformly random between \[-1,1\].

# Generating permutations

Prior to running permutations tests, a separate routine needs to be called to
pre-generate permutations to be used for calculations. This is done to reduce
complexity in the actual calculations program. Simply call the python program
*createPerm.py* to generate permutations:

                python createPerm.py -nPerm <nP> -nA <X> -nB <Y> -o <output permutations txt file>
                
Where *nP* is the number of permutations (e.g 5000), *X* and *Y* are number of
subjects in group A and group B respectively.

Note that the output is the indices of one of the groups (group A) for each
permutation, and NOT one-hot labelling of the groups. This is done to save
space, and permutation calculation routines expects the input to be in this way.

# Running permutation tests

Performing permutation tests can be done simply by calling:

                ./<permutation_test_prog> <file list> <permutation file> <output>

Where *\<permutation_test_prog\>* is either permutationTest_cuda, or
permutationTest_omp in the build folder (assuming cmake is used to compile the
project) for GPU and CPU implementations respectively.

For the CUDA implementation, performing the calculations in parts is supported
when the GPU memory capacity is deemed insufficient to carry out the entire
calculations in one go. This is **NOT** the case for the CPU/OMP implementation.



