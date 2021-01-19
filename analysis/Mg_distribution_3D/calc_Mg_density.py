#!/usr/bin/env python

from math import floor, ceil
import numpy as np
from gridData import Grid

from cafysis.file_io.dcd import DcdFile
from cafysis.file_io.pdb import PdbFile

###############################################################################
# Fixed parameters
NMP_RNA = 1037
CUT_DIST = 15.0
#CUT_DIST2 = 20.0**2
###############################################################################

#print('nMg = {:d}'.format(nMg))


class minmax(object):
    def __init__(self,):
        self.min = [999999.9]*3
        self.max = [-999999.9]*3

    def update(self, xyz):
        for i in range(3):
            if xyz[i] < self.min[i]:
                self.min[i] = xyz[i]
            if xyz[i] > self.max[i]:
                self.max[i] = xyz[i]

    def in_or_out(self, xyz, margin):
        for i in range(3):
            if xyz[i] < (self.min[i] - margin):
                return False
            if xyz[i] > (self.max[i] + margin):
                return False
        return True



if __name__ == '__main__':

    ## H25
    #name = 'h25'
    #cMs  = ['0.0300', '0.0100', '0.0050', '0.0040', '0.0030', '0.0025']
    #resids = list(range(821,830+1)) + list(range(856, 879+1))

    ## h20
    #name = 'h20'
    #cMs  = ['0.0300', '0.0100', '0.0050', '0.0040', '0.0030', '0.0025']
    #resids = list(range(577,586+1)) + list(range(755, 764+1))

    ## CJ
    #name = 'CJ'
    #cMs  = ['0.0300', '0.0100', '0.0050', '0.0040', '0.0030', '0.0025']
    #resids = list(range(569,578+1)) + list(range(817, 821+1))

    # TWJ
    name = 'TWJ'
    cMs  = ['0.0300', '0.0100', '0.0050', '0.0040', '0.0030', '0.0025']
    resids = list(range(586,590+1)) + list(range(649, 655+1)) + list(range(751, 755+1))

    serials = []
    # Use all P, S, B
    for i in sorted(resids):
        if i <= 841:
            res = i - 561
            serials.append(3*(res-1))  # P
            serials.append(3*(res-1)+1)  # S
        elif 848 <= i <= 913:
            res = i - 567
            serials.append(3*(res-1))  # P
            serials.append(3*(res-1)+1)  # S
        else:
            print("Error: resid {:d} should not be there".format(i))
            sys.exit(2)
    
    for cM in cMs:

        filename_dcd = 'cM{:s}.fit{:s}.dcd'.format(cM,name)
        filename_dx  = 'cM{:s}.fit{:s}.dx'.format(cM,name)
    
        nMg = 0
        c = PdbFile('../make_ninfo/16SCD.cM{:s}.cg.ion.pdb'.format(cM),'r').read_all_and_close()[0]
        for r in c.residues:
            if r.atoms[0].name.strip() == 'MG':
                nMg += 1
        
        dcd = DcdFile(filename_dcd)
        dcd.open_to_read()
        dcd.read_header()
        
        data = []
        
        mm = minmax()
        
        ''' Create target volume '''
        mm_target = minmax()
        c = PdbFile('./16SCD.cg.pdb','r').read_all_and_close()[0]
        for i in range(c.num_atom()):
            if i+1 in serials:
                mm_target.update(c.get_atom(i).xyz.get_as_list())
        
        #print('Target volume:')
        #print('    x: {:f} {:f}'.format(mm_target.min[0], mm_target.max[0]))
        #print('    y: {:f} {:f}'.format(mm_target.min[1], mm_target.max[1]))
        #print('    z: {:f} {:f}'.format(mm_target.min[2], mm_target.max[2]))
        
        
        #istep = 0
        while dcd.has_more_data() :
        
            coords_dcd = dcd.read_onestep()
        
            for iMg in range(NMP_RNA, NMP_RNA+nMg):
                xyz = coords_dcd[iMg]
                mm.update(xyz)
                data.append(xyz)
        
            #istep += 1
            #if istep % 1000 == 0:
            #    print(istep)
            #if istep > 2000:
            #    break
        
        #print('{:d} coordinates loaded.'.format(len(data)))
        
        #print('DCD Min-Max:')
        #print('    x: {:f} {:f}'.format(mm.min[0], mm.max[0]))
        #print('    y: {:f} {:f}'.format(mm.min[1], mm.max[1]))
        #print('    z: {:f} {:f}'.format(mm.min[2], mm.max[2]))
        
        #data_target = []
        #for xyz in data:
        #    if mm_target.in_or_out(xyz, CUT_DIST):
        #        data_target.append(xyz)
        
        bins = []
        for i in range(3):
            b = []
            for k in range(floor(mm_target.min[i] - CUT_DIST), ceil(mm_target.max[i] + CUT_DIST)+1):
                b.append(float(k))
                b.append(float(k)+0.5)
            bins.append(b)
        #print (bins)
        
        H, edges = np.histogramdd(np.array(data), bins=bins)
        #print (H.shape, edges[0].size, edges[1].size, edges[2].size)
        #print (H)
        
        g = Grid(H, edges=edges)         
        
        g.export(filename_dx)
