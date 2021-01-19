#!/usr/bin/env python

#from Superimpose_weight import superimpose
from Superimpose_mask import superimpose
from cafysis.file_io.pdb import PdbFile
from cafysis.util_pdb import chains_to_ndarray

NATOM = 1037
chains_ref = PdbFile('../16SCD.cg.pdb','r').read_all_and_close()
chains_que = PdbFile('./cg.pdb','r').read_all_and_close()

#weight = [1.0]*NATOM
weight = [0.0]*NATOM
for i in range(14,43):
    weight[i] = 1.0

d_ref = chains_to_ndarray(chains_ref)
d_que = chains_to_ndarray(chains_que)

rmsd = superimpose(d_ref.T, d_que.T, weight)

print(rmsd)
