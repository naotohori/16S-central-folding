#!/usr/bin/env python

import sys
import math


# Distance where two proton charges have electrostatic energy 2kT
#cut_Mg = 14.5368
#cut_K  =  7.2684
#cut_Cl = 14.5368

#cut =  7.2684
#cut =  7.300000001  # because salt_nt is 0.1 A interval 
cut =  4.40  # 2.1 + 0.8 + 1.5

file_in = open(sys.argv[1],'r')
ntids = []
nPs = []
nMgs = []
nKs  = []
nCls = []
for l in file_in:
    if l[0:6] == '#nstep':
        nstep = int(l.split()[-1])
        print '#',nstep
        continue

    if len(l.strip()) == 0:
        continue
    
    if l[0:1] == '#':
        if len(l.strip()) == 1:
            continue
        ntid = int(l.split()[1])
        ntids.append(ntid)
        nPs.append(0)
        nMgs.append(0)
        nKs.append(0)
        nCls.append(0)
        continue

    lsp = l.split()
    rmin = float(lsp[0])
    rmax = float(lsp[1])
    if rmax <= cut:
        nPs[-1] += int(lsp[6])
        nMgs[-1] += int(lsp[7])
        nKs[-1]  += int(lsp[8])
        nCls[-1] += int(lsp[9])
    

# Use mol/(dm^3) unit
div_factor = float(nstep) * 4/3.0 * math.pi * (cut**3) * 6.02214 / (10**4)

for ntid, nP, nMg, nK, nCl in zip(ntids, nPs, nMgs, nKs, nCls):
    print ntid, nP/div_factor, nMg/div_factor, nK/div_factor, nCl/div_factor
