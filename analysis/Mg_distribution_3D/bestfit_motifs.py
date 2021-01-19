#!/usr/bin/env python

from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.pdb import PdbFile
from Superimpose_mask import superimpose
from cafysis.bestfit_pdb_chains import apply_rot
from numpy import zeros, asarray, float64, matmul
import sys


def dcd_fit(filename_dcd, filename_dcd_out, natom_total, serials, filename_rmsd):
    f_out = open(filename_rmsd,'w')

    # Coord1
    pdb = PdbFile('16SCD.cg.pdb')
    pdb.open_to_read()
    ref_chains = pdb.read_all()
    pdb.close()
    
    #num_atom = 0
    #for chain in ref_chains :
    #    for residue in chain.residues :
    #        num_atom += len(residue.atoms)
        
    ref = zeros((natom_total, 3), dtype=float64, order='C')
    
    i = 0
    for chain in ref_chains :
        for residue in chain.residues:
            for atom in residue.atoms :
                (ref[i][0], ref[i][1], ref[i][2]) = atom.xyz.get_as_tuple()
                i += 1
    
    mask = []
    for i in range(natom_total):
        # the serial ID starts from 1, thus i+1 is the serial ID
        if i+1 in serials:
            mask.append(1)
        else:
            mask.append(0)
    
    dcd = DcdFile(filename_dcd)
    dcd.open_to_read()
    dcd.read_header()
    
    out_dcd = DcdFile(filename_dcd_out)
    out_dcd.open_to_write()
    out_dcd.set_header(dcd.get_header())
    out_dcd.write_header()
    
    #dcd.show_header()
    k = 0
    while dcd.has_more_data() :
        k += 1
        coords_dcd = dcd.read_onestep_np()
    
        rmsd = superimpose(ref.T, coords_dcd.T, mask)
        
        f_out.write('{:8d} {:6.2f}\n'.format(k, rmsd))
    
        out_dcd.write_onestep(coords_dcd)
    
    dcd.close()
    out_dcd.close()

if __name__ == '__main__':
    
    ## H25
    #name = 'h25'
    #cMs  = ['0.0300', '0.0100', '0.0050', '0.0040', '0.0030', '0.0025']
    #nmps = [5594,     4046,    3659,     3581,     3503,     3464] 
    #resids = list(range(821,830+1)) + list(range(856, 879+1))

    ## H20
    #name = 'h20'
    #cMs  = ['0.0300', '0.0100', '0.0050', '0.0040', '0.0030', '0.0025']
    #nmps = [5594,     4046,    3659,     3581,     3503,     3464] 
    #resids = list(range(577,586+1)) + list(range(755, 764+1))

    # CJ
    #name = 'CJ'
    #cMs  = ['0.0300', '0.0100', '0.0050', '0.0040', '0.0030', '0.0025']
    #nmps = [5594,     4046,    3659,     3581,     3503,     3464] 
    #resids = list(range(569,578+1)) + list(range(817, 821+1))

    # TWJ
    name = 'TWJ'
    cMs  = ['0.0300', '0.0100', '0.0050', '0.0040', '0.0030', '0.0025']
    nmps = [5594,     4046,    3659,     3581,     3503,     3464] 
    resids = list(range(586,590+1)) + list(range(649, 655+1)) + list(range(751, 755+1))


    serials = []
    # Use P and S
    for i in sorted(resids):
        if i <= 841:
            res = i - 561
            serials.append(3*(res-1))  # P
            serials.append(3*(res-1)+1)  # S
        elif 848 <= i <= 913:
            res = i - 567
            serials.append(3*(res-1))  # P
            serials.append(3*(res-1)+1)  # S
        else:
            print("Error: resid {:d} should not be there".format(i))
            sys.exit(2)
    #§print(serials)


    for cM, nmp in zip(cMs, nmps):
        dcd_fit('../dcd/cM{:s}.dcd'.format(cM), 'cM{:s}.fit{:s}.dcd'.format(cM,name), nmp, serials, 'cM{:s}.fit{:s}.rmsd'.format(cM, name))

