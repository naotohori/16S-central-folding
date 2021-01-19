#!/usr/bin/env python

from cafysis.elements.pdb import Chain, Residue
from cafysis.file_io.pdb import PdbFile
from cafysis.pdb_aa2cg import aa2cg
from cafysis.bestfit_pdb_chains import fit, apply_rot

# Reference structure 
chains_1j5e_cg = PdbFile('../16SCD.cg.pdb','r').read_all_and_close()


'''
4y4o-pdb-bundle1.pdb
4y4o-pdb-bundle2.pdb  <=== A 562-913 = 16SCD
4y4o-pdb-bundle3.pdb
4y4o-pdb-bundle4.pdb  <=== A 562-913 = 16SCD
'''

DISTANCE_MG = 6.0
def check_distance(chain_ref, atom_query):
    for r in chain_ref.residues:
        for a in r.atoms:
            if atom_query.xyz.distance(a.xyz) <= DISTANCE_MG:
                return True
    return False

for ibundle_16SCD, ibundles_Mg, modelname in ((2, (1,2), '1a'), (4, (3,4), '2a')):

    # Read 16SCD  (Turn on HETATM to read Mg)
    chains = PdbFile('4y4o-pdb-bundle{:d}.pdb'.format(ibundle_16SCD), openmode='r',flg_HETATM=True).read_all_and_close()
    
    c_16SCD = Chain()
    
    for c in chains:
        if c.residues[0].atoms[0].chain_id.strip() != 'A':
            continue
    
        for r in c.residues:
            ''' Some HETATM lines have chain_id 'A' so have to exclude them '''
            if r.atoms[0].chain_id.strip() != 'A':
                continue
    
            if 562 <= r.atoms[0].res_seq <= 913:
                c_16SCD.push_residue(r)
    
    c_Mg = Chain()
    
    # Read Mg
    for ibundle in ibundles_Mg:
        pdb = PdbFile('4y4o-pdb-bundle{:d}.pdb'.format(ibundle), 'r')
        pdb.flg_HETATM = True
        chains = pdb.read_all()
        pdb.close()
        for c in chains:
            for r in c.residues:
                for a in r.atoms:
                    if a.name.strip() == 'MG':
                        if check_distance(c_16SCD, a):
                            r = Residue()
                            r.push_atom(a)
                            c_Mg.push_residue(r)
                            #print(ibundle)
    
    # Only RNA
    pdb_out = PdbFile('4y4o_{:s}_16SCD.original.pdb'.format(modelname), 'w')
    pdb_out.write_all([c_16SCD, ])
    pdb_out.close()

    # RNA and Mg
    pdb_out = PdbFile('4y4o_{:s}_16SCD_Mg.original.pdb'.format(modelname),'w')
    pdb_out.write_all([c_16SCD, c_Mg])
    pdb_out.close()
    
    chains_cg = aa2cg([c_16SCD,])
    rmsd, rot = fit(chains_1j5e_cg, chains_cg)
    apply_rot([c_16SCD, c_Mg], rot)

    # Only RNA
    pdb_out = PdbFile('4y4o_{:s}_16SCD.pdb'.format(modelname), 'w')
    pdb_out.write_all([c_16SCD, ])
    pdb_out.close()

    # RNA and Mg
    pdb_out = PdbFile('4y4o_{:s}_16SCD_Mg.pdb'.format(modelname),'w')
    pdb_out.write_all([c_16SCD, c_Mg])
    pdb_out.close()
