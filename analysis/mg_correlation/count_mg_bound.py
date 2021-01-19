#!/usr/bin/env python

Nmp = 345

nbound = [0]*1037
step = 0
for l in open('test'):
    step += 1
    for i in range(Nmp):
        if int(l[i]) > 0:
            nbound[i] += 1

for i in range(Nmp):
    print ('%i %i %f' % (i+1, nbound[i], nbound[i] / float(step)))

