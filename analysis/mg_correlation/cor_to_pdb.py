#!/usr/bin/env python

import sys

if len(sys.argv) != 3:
    print ('Usage: SCRIPT [input cor file] [output PDB]')
    sys.exit(2)

cor = {}
cor[1] = 0.0
min_cor = 1.0
max_cor = -1.0
for l in open(sys.argv[1],'r'):
    lsp = l.split()
    c = float(lsp[5])
    cor[int(lsp[0])] = c
    if c < min_cor:
        min_cor = c
    if c > max_cor:
        max_cor = c

f_out = open(sys.argv[-1],'w')
for il, l in enumerate( open('../1j5e_SCD.pdb') ):
    resid = int(l[23:26])

    if resid <= 841:
        ic = resid - 561
    else:
        ic = resid - 567

     
    # First two lines, put 1 and -1 to normalize color scale in VMD.
    # Residue 1 (841) does not have Phosphate. (cor data starts from residue 2)
    #if il == 0:
    #    f_out.write('%s%6.2f\n' % (l[0:60], 1.0))
    #elif il == 1:
    #    f_out.write('%s%6.2f\n' % (l[0:60], -1.0))
    #else:
    #    f_out.write('%s%6.2f\n' % (l[0:60], cor[ic]))

    ## Normalize the lower limit to be -max
    #if il == 0:
    #    f_out.write('%s%6.2f\n' % (l[0:60], -max_cor))
    #elif il == 1:
    #    f_out.write('%s%6.2f\n' % (l[0:60], -1.0))
    #else:
    #    f_out.write('%s%6.2f\n' % (l[0:60], cor[ic]))

    # Normalize the upper limit to be -min
    if il == 0:
        f_out.write('%s%6.2f\n' % (l[0:60], -min_cor))
    #elif il == 1:
    #    f_out.write('%s%6.2f\n' % (l[0:60], -1.0))
    else:
        f_out.write('%s%6.2f\n' % (l[0:60], cor[ic]))
