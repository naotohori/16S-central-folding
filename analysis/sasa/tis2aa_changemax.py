#!/usr/bin/env python

import glob
import os
import subprocess


cMs = [0.0, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012, 0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005, 0.01, 0.02, 0.03]

for cM in cMs:

    f_slurm = open('tis2aa_changemax_cM%6.4f.slurm' % (cM,),'w')
    f_slurm.write('#!/bin/bash\n')
    f_slurm.write('#SBATCH -J tis2aa\n')
    f_slurm.write('#SBATCH -n 1\n')
    f_slurm.write('#SBATCH -t 100:00:00\n')
    f_slurm.write('#SBATCH --mem=500m\n')

    cg = glob.glob('cM%6.4f_cg/*pdb' % (cM,))
    aa = glob.glob('cM%6.4f_aa/*pdb' % (cM,))
    
    aabase = []
    for f in aa:
        aabase.append(os.path.basename(f))
    
    for f in cg:
        base = os.path.basename(f)
        if base in aabase:
            pass
        else:
            print base
            f_slurm.write('tis2aa.py --maxpseudo 40 --rmsd 8.0')
            f_slurm.write(' --log cM%6.4f_aa/%s_changemax.log' % (cM, base[:-4]))
            f_slurm.write(' cM%6.4f_cg/%s.pdb' % (cM, base[:-4],))
            f_slurm.write(' cM%6.4f_aa/%s.pdb' % (cM, base[:-4]))
            f_slurm.write('\n')

    f_slurm.close()

    subprocess.call(['sbatch','tis2aa_changemax_cM%6.4f.slurm' % (cM,)])
