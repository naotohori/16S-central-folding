#!/usr/bin/env python

import os
import glob

cMs = [0.0, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012, 0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005, 0.01, 0.02, 0.03]

for cM in cMs:
    pdbs = glob.glob('./cM%6.4f_aamin/*.pdb' % (cM,))

    for p in pdbs:
        base = os.path.basename(p)
    
        f_out = open('./cM%6.4f_aamin_pops/%s' % (cM,base),'w')

        for l in open(p):
            if l[0:4] == 'ATOM':
                seq = l[16:22].strip()
                if len(seq) > 1:
                    seq = seq[0:1]
                f_out.write(l[0:16] + '   %s ' % seq + l[21:])
            else:
                f_out.write(l)
    
    
