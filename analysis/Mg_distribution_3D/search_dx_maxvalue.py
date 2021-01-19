#!/usr/bin/env python

import sys

xmax = 0.0
nlines = 0

for il, l in enumerate(open(sys.argv[1])):
    if il < 12:
        continue
    elif il == 12:
        lsp = l.split()
        nlines = int(lsp[9]) // 3
    elif il < 13+nlines:
        lsp = l.split()
        if len(lsp) != 3:
            print ("Error")
            sys.exit(2)
        lsp = list(map(float, lsp))
        for i in range(3):
            if lsp[i] > xmax:
                xmax = lsp[i]
    else:
        break

print(xmax)
