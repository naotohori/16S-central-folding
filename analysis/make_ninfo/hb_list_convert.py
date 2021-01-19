#!/usr/bin/env python

import sys
#CAN     4  B   19  B  3
#CAN     5  B   18  B  3
#NON     5  B   21  B  1
#NON     5  S   21  B  1  N4 O6 N3 N1 O2 N2

f_out = open('16SCD.hb.list','w')

for l in open('16SCD_16SCD_unprocessed_bonds.dat'):
    lsp = l.strip().split()

    c1 = lsp[0][0:3]

    res1 = int(lsp[1])
    if res1 < 842:
        nt1 = res1 - 561
    else:
        nt1 = res1 - 567
    res2 = int(lsp[2])
    if res2 < 842:
        nt2 = res2 - 561
    else:
        nt2 = res2 - 567

    n_pair = (len(lsp) - 4) / 2

    atoms = []
    types = []
    for ipair in range(n_pair):
        atom1 = lsp[3 + 2*ipair] # column
        if atom1 in ("OP1", "OP2"):
            type1 = 'P'
        elif atom1 in ("O2'", "O4'"):
            type1 = 'S'
        elif atom1 in ("N1", "N2", "N3", "N4", "N6", "N7", "O2", "O4", "O6"):
            type1 = 'B'
        else:
            print ('Error: atom1 type does not exist')
            sys.exit(2)
        atom2 = lsp[3 + 2*ipair + 1] # column
        if atom2 in ("OP1", "OP2"):
            type2 = 'P'
        elif atom2 in ("O2'", "O4'"):
            type2 = 'S'
        elif atom2 in ("N1", "N2", "N3", "N4", "N6", "N7", "O2", "O4", "O6"):
            type2 = 'B'
        else:
            print ('Error: atom2 type does not exist')
            sys.exit(2)

        atoms.append( (atom1, atom2) )
        types.append( (type1, type2) )

    nPP = 0 ; atoms_PP = []
    nPS = 0 ; atoms_PS = []
    nPB = 0 ; atoms_PB = []
    nSP = 0 ; atoms_SP = []
    nSS = 0 ; atoms_SS = []
    nSB = 0 ; atoms_SB = []
    nBP = 0 ; atoms_BP = []
    nBS = 0 ; atoms_BS = []
    nBB = 0 ; atoms_BB = []

    for atm, typ in zip(atoms,types):
        if (typ[0] == 'P' and typ[1] == 'P'):
            nPP += 1
            atoms_PP.append(atm)
        elif (typ[0] == 'P' and typ[1] == 'S'):
            nPS += 1
            atoms_PS.append(atm)
        elif (typ[0] == 'P' and typ[1] == 'B'):
            nPB += 1
            atoms_PB.append(atm)
        elif (typ[0] == 'S' and typ[1] == 'P'):
            nSP += 1
            atoms_SP.append(atm)
        elif (typ[0] == 'S' and typ[1] == 'S'):
            nSS += 1
            atoms_SS.append(atm)
        elif (typ[0] == 'S' and typ[1] == 'B'):
            nSB += 1
            atoms_SB.append(atm)
        elif (typ[0] == 'B' and typ[1] == 'P'):
            nBP += 1
            atoms_BP.append(atm)
        elif (typ[0] == 'B' and typ[1] == 'S'):
            nBS += 1
            atoms_BS.append(atm)
        elif (typ[0] == 'B' and typ[1] == 'B'):
            nBB += 1
            atoms_BB.append(atm)

    if c1 == 'CAN':
        if (nPP > 0 or nPS > 0 or nPB > 0 or
            nSP > 0 or nSS > 0 or nSB > 0 or
            nBP > 0 or nBS > 0 or nBB == 0):
            print ('Error: CAN has to be B-B')
            print (l)
            sys.exit(2)

    if nPP > 0:
        f_out.write('%s %4i P %4i P %i' % (c1, nt1, nt2, nPP))
        for atm in atoms_PP:
            f_out.write(' %s %s' % (atm[0],atm[1]))
        f_out.write('\n')
    if nPS > 0:
        f_out.write('%s %4i P %4i S %i' % (c1, nt1, nt2, nPS))
        for atm in atoms_PS:
            f_out.write(' %s %s' % (atm[0],atm[1]))
        f_out.write('\n')
    if nPB > 0:
        f_out.write('%s %4i P %4i B %i' % (c1, nt1, nt2, nPB))
        for atm in atoms_PB:
            f_out.write(' %s %s' % (atm[0],atm[1]))
        f_out.write('\n')

    if nSP > 0:
        f_out.write('%s %4i S %4i P %i' % (c1, nt1, nt2, nSP))
        for atm in atoms_SP:
            f_out.write(' %s %s' % (atm[0],atm[1]))
        f_out.write('\n')
    if nSS > 0:
        f_out.write('%s %4i S %4i S %i' % (c1, nt1, nt2, nSS))
        for atm in atoms_SS:
            f_out.write(' %s %s' % (atm[0],atm[1]))
        f_out.write('\n')
    if nSB > 0:
        f_out.write('%s %4i S %4i B %i' % (c1, nt1, nt2, nSB))
        for atm in atoms_SB:
            f_out.write(' %s %s' % (atm[0],atm[1]))
        f_out.write('\n')

    if nBP > 0:
        f_out.write('%s %4i B %4i P %i' % (c1, nt1, nt2, nBP))
        for atm in atoms_BP:
            f_out.write(' %s %s' % (atm[0],atm[1]))
        f_out.write('\n')
    if nBS > 0:
        f_out.write('%s %4i B %4i S %i' % (c1, nt1, nt2, nBS))
        for atm in atoms_BS:
            f_out.write(' %s %s' % (atm[0],atm[1]))
        f_out.write('\n')
    if nBB > 0:
        f_out.write('%s %4i B %4i B %i' % (c1, nt1, nt2, nBB))
        for atm in atoms_BB:
            f_out.write(' %s %s' % (atm[0],atm[1]))
        f_out.write('\n')

f_out.close()
