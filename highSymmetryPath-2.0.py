#!/usr/bin/env python
'''
Created on Sep 7, 2016
modifed on Jun 5, 2020 for python3

@author: fu
'''

from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import Poscar
import os

class HighSymmetryPath(object):
    '''
    get the high symmetry path for a given structure
    
    Ags:
        path: path of POSCAR excluding self-POSCAR
    '''

    def __init__(self, path=None):
        if path == None:
            self.path='./'
        else:
            self.path=path
            if not(os.path.exists(self.path)):
                print('Warning! path is wrong.')
            elif not(os.path.isdir(self.path)):
                print("Warning! not a directory.")
                
    def read(self, filename='POSCAR'):
        poscar=Poscar.from_file(self.path+'/'+filename, True, True)
        return poscar.structure
    
    def highSymmetryPath(self, filename='POSCAR', points=40):
        '''
        output the high symmetry path of POSCAR to KPOINTS file
        
        Args:
            filename: POSCAR or CONTCAR
            points: number of kpoints among two high symmetry points
        '''
        
        structure=self.read(filename)
        newStructure=HighSymmKpath(structure)
        kpath=newStructure.kpath['path']
        kpoint=newStructure.kpath['kpoints']
        
        outfile=self.path+'/KPOINTS'
        
        if (os.path.exists(outfile)):
            os.system('rm %s' %outfile)
        with open(outfile, 'a') as out:
            out.write('K-Points\n')
            out.write(' %d\n' %points)
            out.write('Line-mode\n')
            out.write('reciprocal\n')
            
            for i in range(len(kpath)): # number of path
                for k in range(len(kpath[i])): # high symmetry points of a path
                    if (i == 0 and k == 0) or (i == len(kpath) -1 and k == len(kpath[i]) -1):
                        out.write(' %10.6f \t%10.6f \t%10.6f # %s\n' %(kpoint[kpath[i][k]][0], kpoint[kpath[i][k]][1], kpoint[kpath[i][k]][2], kpath[i][k]))
                    else:
                        for j in range(0,2):
                            out.write(' %10.6f \t%10.6f \t%10.6f # %s\n' %(kpoint[kpath[i][k]][0], kpoint[kpath[i][k]][1], kpoint[kpath[i][k]][2], kpath[i][k]))
            
        
# -------------------- test --------------------
path=os.getcwd()
h=HighSymmetryPath(path)
h.highSymmetryPath()
        
