#!/usr/bin/env python

import glob
import os
import subprocess

#f_sh = open('tis2aa.sh','w')
Njob = 200

cgpdbs = glob.glob('./cM0.????_cg/*.pdb')
Ncgpdbs = len(cgpdbs)
Q = Ncgpdbs // Njob
R = Ncgpdbs % Njob

partition = [Q+1 if i < R else Q for i in range(Njob)]

print Ncgpdbs
print partition

N = 0
for i in range(Njob):
    f_slurm = open('tis2aa%03i.slurm' % (i,),'w')
    f_slurm.write('#!/bin/bash\n')
    f_slurm.write('#SBATCH -J tis2aa%03i\n' % (i,))
    f_slurm.write('#SBATCH -n 1\n')
    f_slurm.write('#SBATCH -t 100:00:00\n')
    f_slurm.write('#SBATCH --mem=500m\n')

    for f in cgpdbs[N:N+partition[i]]:
        path = os.path.dirname(f)
        name = os.path.basename(f)[0:-4]
        path_aa = path[0:-3] + '_aa'
        f_slurm.write('tis2aa.py')
        f_slurm.write(' --log %s/%s.log' % (path_aa, name))
        f_slurm.write(' %s' % (f,))
        f_slurm.write(' %s/%s.pdb' % (path_aa, name))
        f_slurm.write('\n')
    f_slurm.close()

    subprocess.call(['sbatch','tis2aa%03i.slurm' % (i,)])

    N = N + partition[i]
