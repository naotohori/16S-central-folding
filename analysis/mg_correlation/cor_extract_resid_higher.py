#!/usr/bin/env python

import sys

if len(sys.argv) != 3:
    print 'Usage: SCRIPT [input cor out] [threshold]'
    sys.exit(2)

thr = float(sys.argv[2])

resids = []
for l in open(sys.argv[1]):
    lsp = l.split()
    res = int(lsp[0])
    if res <= 280:
        resid = res + 561
    else:
        resid = res + 567
    if float(lsp[5]) < thr:
        resids.append(resid)

s = ''
for resid in resids:
    s += ' %i' % (resid,)

print s
