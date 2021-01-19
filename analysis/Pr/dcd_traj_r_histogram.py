#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2015/06/17
@author: Naoto Hori
'''
import sys
import math
import numpy as np
#from cafysis.lib_f2py import py_dcd_r2_histogram_PBC
from cafysis.lib_f2py import py_dcd_r2_histogram
from cafysis.file_io.dcd import DcdFile

nmp = 1037
natom = nmp
npair = natom * (natom-1) / 2
max_r = 350.0
dr = 0.1
nbin = 3500
bin_edges = [x*0.1 for x in range(nbin+1)]
bin_edges_sq = [(x*0.1)**2 for x in range(nbin+1)]



''' Load data of equilibrium Mg = 5 mM to fill trajectory after folded.'''
''' CAUTION: Confirm that the bin_edges are consistent each other '''
#hist_native = np.zeros((nbin,))
#for il, l in enumerate(open('./r_histogram_12_050.out')):
    #hist_native[il] = float(l.split()[2])


hist_all = np.zeros((nbin,))

#concs = [0.0000, 0.0004, 0.0006, 0.0008, 0.0010, 0.0012, 0.0015, 0.0020, 0.0025, 0.0030, 0.0040, 0.0050, 0.0100, 0.0200, 0.0300]
#concs = [0.0000, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.0010, 0.0012, 0.0015, 0.0020, 0.0025, 0.0030, 0.0040, 0.0050, 0.0100, 0.0200, 0.0300]

dcd = DcdFile(sys.argv[1])
dcd.open_to_read()
dcd.read_header()

icount = 0
while dcd.has_more_data() :

    data = dcd.read_onestep_npF()

    #hist = py_dcd_r2_histogram_PBC.dcd_r2_histogram(data[:,:nmp], 350.0, bin_edges_sq)
    hist = py_dcd_r2_histogram.dcd_r2_histogram(data[:,:nmp], bin_edges_sq)
   
    hist_all += hist
    icount += 1

dcd.close()

factor = 1.0 /  float(icount * npair) / dr

f_out = open(sys.argv[2],'w')
for i in range(nbin):
    f_out.write('%f %f %f %f\n' % (bin_edges[i], bin_edges[i+1], hist_all[i] / float(icount), factor * hist_all[i]))
f_out.close()

