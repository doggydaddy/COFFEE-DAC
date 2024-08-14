#!/usr/bin/python

# imports
import numpy as np
import argparse
import random
import sys
from time import perf_counter

parser = argparse.ArgumentParser(
                    prog='generatePermutations.py',
                    description=
                    '''generates indices for N unique permutations
                    of group A of size nA and group B of nB indices.'''
                    )
parser.add_argument('-nPerm', '--numberPermutations', 
                    type=int, 
                    help='number of permutations', 
                    default=5000)
parser.add_argument('-nA', '--numberGroupA', 
                    type=int, 
                    help='number of indices in group A')
parser.add_argument('-nB', '--numberGroupB', 
                    type=int, 
                    help='number of indices in group B')
parser.add_argument('-o', '--outputfile',
                    help='output filepath')

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

nA = int(args.numberGroupA)
nB = int(args.numberGroupB)
nrp = int(args.numberPermutations)
outputfile = args.outputfile

print("creating", str(nrp), "permutations")
print(str(nA), "in one group")
print(str(nB), "in the other group")

def genMfromN(nA, nB, onehot=False):
    '''
    generate M (50) random numbers between [0,N-1] without replacement.
    
    can output indices of one group (let's say the first group), sorted, 
    or the one-hot indices (not sorted, for obvious reasons).
    '''

    N = nA + nB
    M = nA
    i = 0 
    rnd_indices = np.zeros(nA) 
    while i < M:
        a_rnd_num = random.randrange(0, N-1)
        skip_flag = 0
        for j in range(i):
            if rnd_indices[j] == a_rnd_num:
                # this is a duplicate, regen number
                skip_flag = 1
                continue

        if skip_flag == 0:
            rnd_indices[i] = a_rnd_num
            i += 1

    if onehot:
        onehot_output = np.zeros(N, type=np.uint16)
        for k in range(M):
            onehot_output[ int(rnd_indices[k]) ] = 1
        return onehot_output
    else:
        return sorted(rnd_indices)

def genPermutations(nA, nB, nperm):
    '''
    generate nperm permutations with 
    M=nA indices amongst N=nA+nB total indices.
    
    outputs the indices only, 
    and not one-hot (to save space)

    note: There are a lot better ways to do this using standard python
    libraries, but this implementation is sufficiently fast and robust that can
    be directly translated to C/C++ libraries.
    '''

    M = nA
    N = nA + nB
    output = np.zeros((nperm, M))
    print(output.shape)
    p = 0
    while p < nperm:
        redo_flag = 0
        if p == 0: # first iteration

            a_perm = genMfromN(nA, nB)
            output[p, :] = a_perm

        else: # p > 0

            a_perm = genMfromN(nA, nB)

            for pp in range(p):
                matches = 0
                for i in range(M):
                    if output[pp, i] != a_perm[i]:
                        continue
                    else:
                        matches += 1
                if np.sum(matches) == M: # the indices are identical to pp, generate new permutation
                    redo_flag = 1
                    continue 

        if redo_flag != 1: 
            output[p, :] = a_perm
            p += 1

    return output

def combination(n, k):
    '''
    returns n choose k (combinatorics)
    
    computes the combination by formula, 
    which is probably will not translate well to C/C++
    '''
    return int( np.math.factorial(n) / (np.math.factorial(k)*np.math.factorial(n-k) ) )

def choose(n, k):
    '''
    returns n choose k (combinatorics)

    computes the combination by recursion.
    
    this is probably the fastest and most efficient method 
    to do this that can be directly translated to C/C++
    '''
    if k == 0:
        return 1
    else:
        return int((n*choose(n-1, k-1))/k)

print("generating", str(nrp), "permutations")
start = perf_counter()
generated_permutations = genPermutations(nA, nB, nrp)
end = perf_counter()
print("took", str(end-start), "seconds")

print("saving ...")
np.savetxt(outputfile, generated_permutations, fmt='% 4d')
print("all done!")