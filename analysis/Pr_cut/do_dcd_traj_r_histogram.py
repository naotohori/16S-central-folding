#!/usr/bin/env python

import subprocess

concs = [0.0000, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.0010, 0.0012, 0.0015, 0.0020, 0.0025, 0.0030, 0.0040, 0.0050, 0.0100, 0.0200, 0.0300]

for conc in concs:

    slurmfile = 'Pr_cut%6.4f.slurm' % (conc,)
    f_slurm = open(slurmfile, 'w')
    f_slurm.write('#!/bin/bash\n')
    f_slurm.write('#SBATCH -J Pr_cut%6.4f\n' % (conc,))
    f_slurm.write('#SBATCH -n 1\n')
    f_slurm.write('#SBATCH -t 100:00:00\n')
    f_slurm.write('#SBATCH --mem=500m\n')
 
    f_slurm.write('\n')
    f_slurm.write('mkdir /tmp/NH_16SCD_Pr_cut_%6.4f\n' % (conc,))
    f_slurm.write('cp ../dcd/cM%6.4f.body.dcd /tmp/NH_16SCD_Pr_cut_%6.4f\n' % (conc,conc,))

    f_slurm.write('./dcd_traj_r_histogram.py /tmp/NH_16SCD_Pr_cut_%6.4f/cM%6.4f.body.dcd ./cM%6.4f.Pr\n' % (conc,conc,conc))

    f_slurm.write('rm -rf /tmp/NH_16SCD_Pr_cut_%6.4f\n' % (conc,))

    f_slurm.close()

    subprocess.call(["sbatch", slurmfile])
