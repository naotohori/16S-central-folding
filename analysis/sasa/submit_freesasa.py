#!/usr/bin/env python

import glob
import os
import subprocess

#f_sh = open('tis2aa.sh','w')
Njob = 200

pdbs = glob.glob('./cM0.????_aamin_pops/*.pdb')

Npdbs = len(pdbs)
Q = Npdbs // Njob
R = Npdbs % Njob

partition = [Q+1 if i < R else Q for i in range(Njob)]

print Npdbs
print partition

N = 0
for i in range(Njob):
    f_slurm = open('frss%03i.slurm' % (i,),'w')
    f_slurm.write('#!/bin/bash\n')
    f_slurm.write('#SBATCH -J frss%03i\n' % (i,))
    f_slurm.write('#SBATCH -n 1\n')
    f_slurm.write('#SBATCH -t 100:00:00\n')
    f_slurm.write('#SBATCH --mem=500m\n')

    for f in pdbs[N:N+partition[i]]:
        path = os.path.dirname(f)
        name = os.path.basename(f)[0:-4]
        f_slurm.write('freesasa')
        f_slurm.write(' -f pdb %s' % (f,))
        f_slurm.write(' -o %s_frss/%s.pdb' % (path[:-11], name))
        f_slurm.write(' -e %s_frss/%s.err' % (path[:-11], name))
        f_slurm.write('\n')
    f_slurm.close()

    subprocess.call(['sbatch','frss%03i.slurm' % (i,)])

    N = N + partition[i]
