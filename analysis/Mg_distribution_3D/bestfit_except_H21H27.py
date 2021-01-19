#!/usr/bin/env python

from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.pdb import PdbFile
from Superimpose_mask import superimpose
from numpy import zeros, asarray, float64
import sys

id_begin = []
id_end = []
for iarg in range(4, len(sys.argv)-1, 2) :
    id_begin.append(int(sys.argv[iarg]))
    id_end.append(int(sys.argv[iarg+1]))

NATOM_TOTAL = 3659
#NATOM_TOTAL = 5594

# Coord1
pdb = PdbFile('16SCD.cg.pdb')
pdb.open_to_read()
ref_chains = pdb.read_all()
pdb.close()

num_atom = 0
for chain in ref_chains :
    for residue in chain.residues :
        num_atom += len(residue.atoms)
    
#ref = zeros((3, NATOM_TOTAL), dtype=float64, order='F')
ref = zeros((NATOM_TOTAL, 3), dtype=float64, order='C')

i = 0
for chain in ref_chains :
    for residue in chain.residues:
        for atom in residue.atoms :
            #(ref[0][i], ref[1][i], ref[2][i]) = atom.xyz.get_as_tuple()
            (ref[i][0], ref[i][1], ref[i][2]) = atom.xyz.get_as_tuple()
            i += 1

# 1 - 77 and 270 - 950
mask = []
for i in range(NATOM_TOTAL):
    if i in range(0,77):
        mask.append(1)
    elif i in range(270-1,950):
        mask.append(1)
    else:
        mask.append(0)

dcd = DcdFile('dcd/cM0.0050.dcd')
dcd.open_to_read()
dcd.read_header()

out_dcd = DcdFile('cM0.0050.fit.dcd')
out_dcd.open_to_write()
out_dcd.set_header(dcd.get_header())
out_dcd.write_header()

out_rmsd = open('cM0.0050.fit.rmsd', 'w')

#dcd.show_header()
k = 0
while dcd.has_more_data() :
    k += 1
    coords_dcd = dcd.read_onestep_np()
    
    rmsd = superimpose(ref.T, coords_dcd.T, mask)
    
    out_rmsd.write('{:f}\n'.format(rmsd))
    out_dcd.write_onestep(coords_dcd)
        

dcd.close()
out_rmsd.close()
out_dcd.close()
