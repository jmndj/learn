#!/usr/bin/env python
#coding=utf-8
'''
Created on Jun. 5, 2018

add SOC [Jul. 6, 2018]

add splitting calculations [Jan. 6, 2019]

fix bug in finding band sides [Feb. 17, 2019]
add spin-related codes [Aug. 16, 2019]

support y-orbit [March, 6, 2022]
plot band by EIGENVAL file

@author: fu
'''
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import math
import collections
from copy import deepcopy
from progressbar import ProgressBar, Percentage, Bar, ETA
from matplotlib.font_manager import FontProperties
import matplotlib.gridspec as gridspec


class PlotBand():
    def __init__(self):
        import argparse
        import yaml
        
        parser=argparse.ArgumentParser(description="plot Band and DOS")
        parser.add_argument('--isSemi', 
                            action='store_true', 
                            help='swtich to semiconductor material (Default: metal). If true, need to shit the Fermi level to the VBM when plotting band.')
        parser.add_argument('--isDecomposed', 
                            action='store_false', 
                            help='swtich from decompose atomic orbits to total orbits')
        
# =============================================================================
#         parser.add_argument('-p', '--path', 
#                             type=str,
#                             nargs='+',
#                             default=os.getcwd(),
#                             help='path of electronic structure data (Default: current directory)')
# =============================================================================
        parser.add_argument('-c', '--calculation', 
                            choices=['normal', 'spin', 'soc', 'spin+soc', 
                                     'hse', 'hse+spin', 'hse+soc', 'hse+spin+soc'], 
                            default='normal',
                            help="type of calculation (Default: 'normal')")
        

        parser.add_argument('-pb', '--path_of_band', type=yaml.load, default=[os.getcwd()], help='path of electronic structure data (Default: current directory)')
        parser.add_argument('-pd', '--path_of_dos', type=str, default=os.getcwd(), help='path of dos (Default: current directory)')
        # parser.add_argument('-e', '--erange', type=yaml.load, default=[-4, 4], help='range of energies')
        parser.add_argument('-e', '--erange', type=float, nargs='+', default=[-4, 4], help='range of energies')
        
        parser.add_argument('-B', '--Band', action='store_true', help='plot band structure (default: False)')
        parser.add_argument('-D', '--PDOS', action='store_true', help='plot PDOS (default: False)')
        parser.add_argument('-BD', '--Band_PDOS', action='store_true', help='plot band structure and its DOS (default: False)')
        
        
        args=parser.parse_args()
        self.isSemi=args.isSemi
        self.isDecomposed=args.isDecomposed
        # self.path=args.path
        self.calculation=args.calculation
        
        self.path_of_band=args.path_of_band
        self.path_of_dos=args.path_of_dos
        self.erange=args.erange
        
        self.Band=args.Band
        self.PDOS=args.PDOS
        self.Band_PDOS=args.Band_PDOS
        print(args)
        if self.Band_PDOS:
            if self.calculation in ['normal']:
                self.path_of_dos='{}/../scf'.format(self.path_of_band[0])
            elif self.calculation in ['soc']:
                self.path_of_dos='{}/../soc'.format(self.path_of_band[0])
            elif args.path_of_dos == os.getcwd():
                raise ValueError('input the path of dos')
        
        
        self.btype={'normal': 0,
                    'spin': 1,
                    'soc': 0,
                    'spin+soc': 0,
                    'hse':0,
                    'hse+spin': 1,
                    'hse+soc': 0,
                    'hse+spin+soc': 0}
        
        self.atype={'normal': 0,
                    'spin': 1,
                    'soc': 2,
                    'spin+soc': 2,
                    'hse':0,
                    'hse+spin': 1,
                    'hse+soc': 2,
                    'hse+spin+soc': 2}
        
        self.ktype={'normal': 0,
               'spin': 0,
               'soc': 0,
               'spin+soc': 0,
               'hse':1,
               'hse+spin': 1,
               'hse+soc': 1,
               'hse+spin+soc': 1}
        
    
    def _getSkipforHSE(self, path, filename='KPOINTS'):
        """
        get number of skipping for plot HSE's band
        """
        infile=open('{}/{}'.format(path, filename))
        
        for i in range(0, 3): infile.readline()
        
        str0=infile.readline()
        counter=0
        while str0:
            w0=int(str0.split()[3])
            if w0 == 0: return counter
            
            str0=infile.readline()
            counter += 1
        return 0
    
    def getReducedBand(self, path, filename='EIGENVAL', skip=0):
        """
        read band data from EIGENVAL file
        
        filename: EIGENVAL file
        skip: skip give kpoints. Note that the counter is from 1.
        
        an example of EIGENVAL's format:
                    electrons(NE)    kpoints(NK)    valence bands(NVB)    Occupancy(noc)
        PBE:          56               30              28                   1   <- btype 0
        PBE+spin:     56               30              28                   2   <- btype 1
        PBE+soc:      56               30              56                   1   <- btype 0
        PBE+spin+soc: 56               30              56                   1   <- btype 0
        
        Return:
            {'kpoints': [nkpoints, 3],
             'bands': [nkpoints, nbands] (non-spin)
                      [nkpoints, nbands, 2] (spin)}
        """
        btype=self.btype
        
        if 'hse' in self.calculation: skip=self._getSkipforHSE(path=path)
        
        infile=open('{}/{}'.format(path, filename))
        for i in range(5): infile.readline()
        
        [nelectrons, nkpoints, nbands]=[int(s0) for s0 in infile.readline().split()]
        print('nelectrons: {}\nnkpoints: {}\nnbands: {}'.format(nelectrons, nkpoints, nbands))
        print('skip: {}'.format(skip))
        
        kpoints=np.empty((nkpoints-skip, 3))
        # if self.calculation in ['normal', 'soc', 'spin+soc', 'hse', 'hse+soc', 'hse+spin+soc']:
        if btype[self.calculation] == 0:
            bands=np.empty((nkpoints-skip, nbands))
        # elif self.calculation in ['spin', 'hse+spin']:
        elif btype[self.calculation] == 1:
            bands=np.empty((nkpoints-skip, nbands, 2))
            
        counter=1 # for kpoints
        for i in range(0, nkpoints):
            infile.readline() # blank line
            if counter > skip:
                kpoints[i-skip][:]=[float(s0) for s0 in infile.readline().split()][:3] # [x, y, z, weight]
            else:
                infile.readline()
                
            bands_of_kpoint0=[]
            for j in range(0, nbands):
                data0=[]
                if counter > skip:
                    # if self.calculation in ['normal', 'soc', 'spin+soc']:
                    if btype[self.calculation] == 0:
                        bands[i-skip][j]=[float(s0) for s0 in infile.readline().split()][1]
                    # elif self.calculation in ['spin', 'hse+spin']:
                    elif btype[self.calculation] == 1:
                        bands[i-skip][j][:]=[float(s0) for s0 in infile.readline().split()][1:3]
                else:
                    infile.readline()
            counter += 1
        
        return {'kpoints': np.array(kpoints), 'bands': np.array(bands)}
    
    def getFullBand(self, path, filename='PROCAR', skip=0):
        """
        read full band datat containg atmoic orbits information from PROCAR file
        
        filename: EIGENVAL file
        skip: skip give kpoints. Note that the counter is from 1.
        
        an example of PROCAR's format:
                    kpoints(NK)    bands(NK)    ions(N)    nstates(PDOS)
        PBE:          30             36           8          x   <- atype 0
        PBE+spin:     30+30          54           8          x   <- atype 1
        PBE+soc:      30             81           8+8+8+8    x   <- atype 2
        PBE+spin+soc: 30             81           8+8+8+8    x   <- atype 2
        
        Return:
            {'kpoints': [nkpoints, 3],
             'bands': [nkpoints, nbands] (non-spin) 
                      [nkpoints, nbands, 2] (spin),
                      
             'atoms': [nkpoints, nbands, natoms, norbits] (normal) 
                      [nkpoints, nbands, natoms, norbits, 2] (spin). Here, 2 -> [up, down] for spin
                      [nkpoints, nbands, natoms, norbits, 4] (soc or spin+soc)}. Here, 4 -> [m_total, m_x, m_y, m_z] for magnettism
        """
        atype=self.atype
        
        if 'hse' in self.calculation: skip=self._getSkipforHSE()
        
        infile=open('{}/{}'.format(path, filename))
        infile.readline()
        
        [nkpoints, nbands, natoms]=[int(s0) for s0 in infile.readline().split() if s0.isdigit()]
        norbits=len(os.popen('head -8 {}/{} | tail -1'.format(self.path, filename)).readline().split())-1
        print('nkpoints: {}\nnbands: {}\nnatoms: {}\nnorbits: {}'.format(nkpoints, nbands, natoms, norbits))
        print('skip: {}'.format(skip))
        
        kpoints=np.empty((nkpoints-skip, 3))
        
        repeat_of_kpoints=1
        repeat_of_atoms=1
        
        if atype[self.calculation] == 0:
            bands=np.empty((nkpoints-skip, nbands))
            atoms=np.empty((nkpoints-skip, nbands, natoms+1, norbits))
        elif atype[self.calculation] == 2:
            bands=np.empty((nkpoints-skip, nbands))
            atoms=np.empty((nkpoints-skip, nbands, natoms+1, norbits, 4))
            repeat_of_atoms=4
        elif atype[self.calculation] == 1:
            bands=np.empty((nkpoints-skip, nbands, 2))
            atoms=np.zeros((nkpoints-skip, nbands, natoms+1, norbits, 2))
            repeat_of_kpoints=2
        
        for rk in range(0, repeat_of_kpoints): # 1: non-spin or spin+soc; 2: spin
            counter=1 # for kpoints
            for i in range(0, nkpoints):
                infile.readline() # blank line
                if counter > skip:
                    str0=infile.readline().split()
                    while str0 == []:
                        str0=infile.readline().split()
                    kpoints[i-skip][:]=[float(s0) for s0 in str0[3:6]] # [x, y, z]
                else:
                    infile.readline()
                    
                infile.readline()
                for j in range(0, nbands):
                    if counter > skip:
                        str0=infile.readline().split()
                        while str0 == []:
                            str0=infile.readline().split()
                        b0=float(str0[4])
                        if repeat_of_kpoints == 1:
                            bands[i-skip][j]=b0
                        else:
                            bands[i-skip][j][rk]=b0
                    else:
                        infile.readline()
                
                    infile.readline()
                    infile.readline()
                    for ra in range(0, repeat_of_atoms): # 1: non-soc; 4: soc [mtotal, mx, my, mz]
                        for k in range(0, atoms.shape[2]): # natoms
                            if counter > skip:
                                str0=infile.readline().split()
                                while str0 == []:
                                    str0=infile.readline().split()
                                data0=[float(s0) for s0 in str0[1:]]
                                if (repeat_of_kpoints == 1) and (repeat_of_atoms == 1):
                                    atoms[i-skip][j][k][:]=data0
                                elif (repeat_of_kpoints > 1) and (repeat_of_atoms == 1):
                                    atoms[i-skip, j, k, :, rk]=data0
                                elif (repeat_of_kpoints == 1) and (repeat_of_atoms > 1):
                                    # atoms[i][j][k][:][ra]=data0[n] # error
                                    atoms[i-skip, j, k, :, ra]=data0
                                else:
                                    raise ValueError('repeat_of_kpoints: {} repeat_of_atoms: {}'.format(repeat_of_kpoints, repeat_of_atoms))
                            else:
                                infile.readline()
                                
                    infile.readline()
                counter += 1
            infile.readline()                
        return {'kpoints': np.array(kpoints), 'bands': np.array(bands), 'atoms': atoms}
    
    def getDOS(self, path, filename='DOSCAR'):
        """
        an example of DOSCAR's format:
                    NEDOS lines                                   PDOS
        PBE:          DOS intg_DOS                                  s px py pz ...                                 <- btype 0 for 1st column
        PBE+spin:     DOS(up) DOS(dn) intg_DOS(up)  intg_DOS(dn)    s(up) s(dn) px(up) px(dn) ...                  <- btype 1
        PBE+soc:      DOS intg_DOS                                  s(total) s(mx) s(my) s(mz) p(total) p(mx) ...  <- btype 0
        PBE+spin+soc: DOS intg_DOS                                  s(total) s(mx) s(my) s(mz) p(total) p(mx) ...  <- btype 0
        """
        btype=self.btype
        
        infile=open('{}/{}'.format(path, filename))
        
        for i in range(0, 5): infile.readline()
        
        str0=infile.readline().split()
        npoints=int(str0[2])
        Efermi=float(str0[3])
        
        # read total DOS
        energies=np.empty(npoints)
        if btype[self.calculation] == 0: # non-spin
            totalDOS=np.empty((npoints))
            cumulativeDOS=np.empty((npoints))
        elif btype[self.calculation] == 1: # spin
            totalDOS=np.empty((npoints, 2))
            cumulativeDOS=np.empty((npoints, 2))
            
        for i in range(0, npoints):
            data0=[float(s0) for s0 in infile.readline().split()]
            e0=data0[0]-Efermi
            
            energies[i]=e0
            if btype[self.calculation] == 0: # non-spin
                totalDOS[i]=data0[1]
                cumulativeDOS[i]=data0[2]
            elif btype[self.calculation] == 1: # spin
                totalDOS[i][0]=data0[1] # up
                totalDOS[i][1]=data0[2] # down
                cumulativeDOS[i][0]=data0[3] # up
                cumulativeDOS[i][1]=data0[4] # down
        
        symbols, elementNumbers, totalAtomNumbers=self.getElementinfo(path=path)
        
        atoms=[] # PDOS per atom
        for i in range(0, totalAtomNumbers):
            infile.readline()
            data0=[]
            for j in range(0, npoints):
                data0.append([float(s0) for s0 in infile.readline().split()][1:])
            atoms.append(data0)
        atoms=np.array(atoms)
        elements={} # PDOS per element
        for s0 in symbols:
            index0=self.getIndexOfAtoms(path=path, symbol_of_element=s0)
            elements[s0]=atoms[index0[0]:index0[1], :, :].sum(axis=0)
        
        return {'energies': energies, 'totalDOS': totalDOS, 'cumulativeDOS': cumulativeDOS, 'atoms': atoms, 'elements': elements}

    def getLatticeVector(self, path, filename='POSCAR'):
        """
        get the lattice vector and reciprocal lattice vector.
        """
        infile=open('{}/{}'.format(path, filename))
        infile.readline() # comment line

        scale=float(infile.readline())
        
        a=[] # lattice vector
        for i in range(0, 3):
            a.append([float(s0)*scale for s0 in infile.readline().split()])
        a=np.array(a)

        volume=a[0][0]*a[1][1]*a[2][2]+a[0][1]*a[1][2]*a[2][0]+ \
               a[0][2]*a[1][0]*a[2][1]-a[0][0]*a[1][2]*a[2][1]- \
               a[0][1]*a[1][0]*a[2][2]-a[0][2]*a[1][1]*a[2][0]
        
        b=np.zeros((3,3)) # reciprocal lattice vector
        for i in range(0,3):
            if i == 0:
                j=1
                k=2
            elif i == 1:
                j=2
                k=0
            else:
                j=0
                k=1
            c=np.zeros(3)
            c[0]=a[j][1]*a[k][2]-a[j][2]*a[k][1]
            c[1]=a[j][2]*a[k][0]-a[j][0]*a[k][2]
            c[2]=a[j][0]*a[k][1]-a[j][1]*a[k][0]
            for j in range(0, 3):
                b[i][j]=2*math.pi*c[j]/volume
       
        return a, b, volume

    def getEfermi(self, path, filename='OUTCAR'):
        """
        get Fermi energy level from OUTCAR.
        """
        return float(os.popen("grep 'E-fermi' {}/{}".format(path, filename)).readline().split()[2])
    
    def getKpointsNum(self, path, filename='KPOINTS'):
        """
        get number of kpoints from KPOINTS file.
        
        Note that: remove the blank line at end of the file if exist.
        """
        ktype=self.ktype
        
        infile='{}/{}'.format(path, filename)
        
        kpoints=None
        
        num_per_path=int(os.popen('head -2 {} | tail -1'.format(infile)).readline())
    
        if ktype[self.calculation] == 0:
            num_of_kpoints=int(os.popen('wc -l {}'.format(infile)).readline().split()[0])-4
            kpoints=num_per_path*num_of_kpoints/2
        elif ktype[self.calculation] == 1:
            num_of_kpoints=int(os.popen('wc -l {}'.format(infile)).readline().split()[0])-3
            kpoints=num_of_kpoints-self._getSkipforHSE()

        return kpoints
        
    def getBandSides(self, bands, Efermi, tolerance=1e-3):
        """
        get band side (VBM and CBM).
        
        Arguments:
            bands: bands [nkpoints, nbands].
            Efermi: Fermi energy level.
            tolerance (default=1e-3): tolerance for Efermi at downside.
            
        Return:
            
        """
        upside=[]
        downside=[]
        
        for i in range(0, bands.shape[0]): # kpoint
            for j in range(0, bands.shape[1]): # band
                if bands[i][j]-Efermi < tolerance: # downside
                    if (downside == []) or ((abs(bands[i][j]-downside[-1]) < tolerance) and (j > downside[1])) or (bands[i][j] > downside[-1]) : 
                        downside=[i, j, bands[i][j]]
                if bands[i][j] > Efermi: # upside
                    if (upside == []) or ((abs(bands[i][j]-downside[-1]) < tolerance) and (j < downside[1]))or (bands[i][j] < upside[-1]): 
                        upside=[i, j, bands[i][j]]

        return upside, downside
    
    def getDirectBandgap(self, bands, Efermi, tolerance=1e-3):
        """
        get the direct bandgap.
        
        Arguments:
            bands: bands [nkpoints, nbands].
            Efermi: Fermi energy level.
            
        Return:
            [nkpoint, nvband, nncband, Eg]
        """
        bandgap=[]
        for i in range(0, bands.shape[0]): # kpoint
            upside=[]
            downside=[]
            for j in range(0, bands.shape[1]): # band
                if bands[i][j]-Efermi < tolerance: # downside
                    if (downside == []) or ((abs(bands[i][j]-downside[-1]) < tolerance) and (j > downside[1])) or (bands[i][j] > downside[-1]) : 
                        downside=[i, j, bands[i][j]]
                if bands[i][j] > Efermi: # upside
                    if (upside == []) or ((abs(bands[i][j]-downside[-1]) < tolerance) and (j < downside[1]))or (bands[i][j] < upside[-1]): 
                        upside=[i, j, bands[i][j]]

            deltaD=upside[2]-downside[2]
            if (bandgap == []) or (bandgap[-1] > deltaD):
                bandgap=[upside[0], upside[1], downside[1], deltaD]

        DirectBandgap=bandgap
        return DirectBandgap
    
    def calcKpointsDistance(self, kpoints, path):
        """
        calculate the distance of two kpoints
        """
        kpointsDistance=[]
        for i in range(0, kpoints.shape[0]):
            if i == 0:
                kpointsDistance=np.array([0.0])
            else:
                c1=kpoints[i]
                c0=kpoints[i-1]

                a, b, volume=self.getLatticeVector(path=path)
                dx=(c1[0]-c0[0])*b[0][0]+(c1[1]-c0[1])*b[1][0]+(c1[2]-c0[2])*b[2][0]
                dy=(c1[0]-c0[0])*b[0][1]+(c1[1]-c0[1])*b[1][1]+(c1[2]-c0[2])*b[2][1]
                dz=(c1[0]-c0[0])*b[0][2]+(c1[1]-c0[1])*b[1][2]+(c1[2]-c0[2])*b[2][2]
                delta=math.sqrt(dx*dx+dy*dy+dz*dz)

                temp=kpointsDistance[i-1]+delta
                kpointsDistance=np.hstack([kpointsDistance, temp])
        
        return kpointsDistance

    def getHighSymmetryPoints(self, path, filename='KPOINTS'):
        """
        get high symmetry points from KPOINTS
        """
        ktype=self.ktype
        
        infile=open('{}/{}'.format(path, filename))
        
        # kpoints and kpoints' symbol
        hspoints=[]
        hspointsSymbols=[]
        
        str0=infile.readline()
        isNeeded=False
        while str0:
            if str0.lower().startswith('r'): 
                isNeeded=True
                str0=infile.readline()
            if isNeeded:
                data0=str0.split()
                if ((ktype[self.calculation] == 0) and (len(data0) == 5)) or ((ktype[self.calculation] == 1) and (len(data0) == 6)):
                    k0=[float(s0) for s0 in data0[:3]]
                    symbol0=data0[-1]
                    if (symbol0.lower() == 'gamma') or (symbol0.lower() == '\gamma'): symbol0='$\Gamma$'
                    
                    hspoints.append(k0)
                    hspointsSymbols.append(symbol0)
                    
            str0=infile.readline()

        hspoints=np.array(hspoints)
        hspointsSymbols=np.array(hspointsSymbols)

        # high symmerty point distance
        hspointsDistances=self.calcKpointsDistance(kpoints=hspoints, path=path)

        return hspointsDistances, hspointsSymbols

    def getElementinfo(self, path, filename='POSCAR'):
        """
        get element information from POSCAR
        """
        
        elements=os.popen('sed -n 6p {}/{}'.format(path, filename)).readline().split()
        string=os.popen('sed -n 7p {}/{}'.format(path, filename)).readline()
        elementNumbers=[int(s0) for s0 in string.split()]
        
        totalAtomNumbers=0
        for i in range(0, len(elementNumbers)):
            totalAtomNumbers += elementNumbers[i]

        return elements, elementNumbers, totalAtomNumbers

    def getIndexOfAtoms(self, path, symbol_of_element):
        """
        get the index of atoms for given element
        
        Arguments:
            element: symbol of element. e.g., 'Li'
            
        Return:
            [start, end]
        """
        elements, elementNumbers, totalAtomNumbers=self.getElementinfo(path=path)
        index_of_atoms=[]
        for i in range(0, len(elements)):
            if elements[i] == symbol_of_element:
                if i == 0:
                    index_of_atoms=np.array([0, elementNumbers[i]])
                else:
                    temp=0
                    for j in range(0, i):
                        temp += elementNumbers[j]
                    index_of_atoms=np.array([temp, elementNumbers[i]+temp])
        return index_of_atoms
    
    def getSatesOfElement(self, symbol_of_element, atoms, nkpoint, nband, norbit, path, **kwargs):
        """
        get the states of given element at (nband, nkpoint, norbit)
        
        Note that: ions don't include the degree of freedom of magnetization
        
        Arguments:
            atoms: projected density of states for all atoms, e.g.,
                    [nkpoints, nbands, natoms, norbits] (normal) 
                    [nkpoints, nbands, natoms, norbits, 2] (spin). Here, 2 -> [up, down] for spin
                    [nkpoints, nbands, natoms, norbits, 4] (soc or spin+soc)}. Here, 4 -> [m_total, m_x, m_y, m_z] for magnettism
            
            nkpoint: number of kpoint
            nband: number of band
            element: symbol of element
            norbit: number of orbit
            
            kwargs:
                for spin:
                    nspin:
                
                for soc:
                    nm:
        """
        atype=self.atype
        
        index_of_atoms=self.getIndexOfAtoms(path=path,
                                            symbol_of_element=symbol_of_element)
        
        state=0
        for i in range(index_of_atoms[0], index_of_atoms[1]):
            if atype[self.calculation] == 0:
                state += atoms[nkpoint][nband][i][norbit] # total states
            elif atype[self.calculation] == 1:
                if 'nspin' not in kwargs: raise ValueError('Need nspin')
                state += atoms[nkpoint][nband][i][norbit][kwargs['nspin']]
            elif atype[self.calculation] == 2:
                if 'nm' not in kwargs: raise ValueError('Need nm')
                state += atoms[nkpoint][nband][i][norbit][kwargs['nm']]
                
        return state
    
    def showBandInfo(self, bands, Efermi, tolerance4bandside):
        """
        """
        upside, downside=self.getBandSides(bands=bands, Efermi=Efermi, tolerance=tolerance4bandside)
        DirectBandgap=self.getDirectBandgap(bands=bands, Efermi=Efermi, tolerance=tolerance4bandside)

        print('\n\n++++++++++++++++++++ Indirect Band ++++++++++++++++++++')
        print('  CBM(nk,nb,E): [{:d}, {:d}, {:.4f}]'.format(upside[0], upside[1], upside[2]))
        print('  VBM(nk,nb,E): [{:d}, {:d}, {:.4f}]'.format(downside[0], downside[1], downside[2]))
        print('  Bandgap_i: {:4.4f}'.format(upside[2]-downside[2]))
        print('+++++++++++++++++++++ Direct Band +++++++++++++++++++++')
        
        print('  Bandgap_d(nk,nvb,ncb,E): [{:d}, {:d}, {:d}, {:.4f}]\n'.format(DirectBandgap[0], DirectBandgap[1], DirectBandgap[2], DirectBandgap[-1]))
        if upside[0] == downside[0]:
            print('>>>>>>>>>>>>>>>>>>> Direct bandgap <<<<<<<<<<<<<<<<<<<<\n\n')
        else:
            print('>>>>>>>>>>>>>>>>>> Indirect bandgap <<<<<<<<<<<<<<<<<<<\n\n')
        return upside, downside, DirectBandgap
    
    def _plotHighSymmetryLines(self, path, plt):
        """
        """
        hspointsDistances=[] 
        hspointsSymbols=[]
        if isinstance(path, str):
            hspointsDistances, hspointsSymbols=self.getHighSymmetryPoints(path=path)
            
        elif isinstance(path, list) or isinstance(path, np.ndarray):
            for p0 in path:
                d0, s0=self.getHighSymmetryPoints(path=p0)
                d0=d0.tolist()
                s0=s0.tolist()
                if hspointsDistances == []:
                    hspointsDistances=d0
                    hspointsSymbols=s0
                else:
                    hspointsDistances += [x+ hspointsDistances[-1] for x in d0]
                    hspointsSymbols += s0
            hspointsDistances=np.array(hspointsDistances)
            hspointsSymbols=np.array(hspointsSymbols)
        for i in range(0, len(hspointsDistances)):
            plt.axvline(x=hspointsDistances[i], color="k")
        
        return hspointsDistances, hspointsSymbols

# =============================================================================
#     def plotBand(self, path_of_structrues, erange, filename, **kwargs):
#         """
#         plot band structure from EIGENVAL file.
#         
#         Arguments:
#             path_of_structures:
#             filename: output filename. e.g. bandfat-Li
#             erange: range of plotted energy. e.g. [-4, 4]
#             
#             kwargs:
#                 isFilledBandGap (default=False): whether to fill the band gap.
#                 
#                 for non-spin:
#                     isFilledBandGap (default: False):
#                     color4FilledBandGap (default: 'b'): the color of filling band gap.
#                 for spin:
#                     isFilledBandGap_up (default: False):
#                     isFilledBandGap_dn (default: False):
#                     color4FilledBandGap_up (default='m'): color of filling band gap.
#                     color4FilledBandGap_dn (default='c'): color of filling band gap.
#                     
#                 tolerance4bandside (default=1e-3):
#                 tolerance4kpoint (defult=1e-6):
#         """
#         btype=self.btype
#         
#         tolerance4bandside=1e-6
#         if 'tolerance4bandside' in kwargs: tolerance4bandside=kwargs['tolerance4bandside']
# # =============================================================================
# #         tolerance4kpoint=1e-6
# #         if 'tolerance4kpoint' in kwargs: tolerance4kpoint=kwargs['tolerance4kpoint']
# # =============================================================================
#         
#         plt.figure(figsize=(8,8)) # size of canvas: 8x8 inch
#         
#         paths=[]
#         if isinstance(path_of_structrues, str): 
#             paths=[path_of_structrues]
#         elif isinstance(path_of_structrues, list) or isinstance(path_of_structrues, np.ndarray):
#             paths=path_of_structrues
#         
#         # prepare data
#         distances_of_all=[]
#         for ip in range(0, len(paths)):
#             path=paths[ip]
#             data=self.getReducedBand(path=path)
#             kpoints=data['kpoints']
#             bands=data['bands']
#         
#             distances=self.calcKpointsDistance(kpoints, path)
#             
#             # process distances
#             if ip == 0:    
#                 distances_of_all=distances.tolist()
#             else:
#                 distances += distances_of_all[-1]
#                 distances_of_all=distances_of_all+distances.tolist()
#                     
#             # plot band
#             Efermi=self.getEfermi(path)
#             if btype[self.calculation] == 0:
#                 upside, downside, DirectBandgap=self.showBandInfo(bands, Efermi, tolerance4bandside=tolerance4bandside)
#                 # move Fermi to VBM
#                 bands -= downside[-1]
#                 for i in range(0, bands.shape[1]): plt.plot(distances, bands[:, i], 'k', lw=1)
#             
#                 # fill bandagap
#                 isFilledBandGap=False
#                 if 'isFilledBandGap' in kwargs: isFilledBandGap=kwargs['isFilledBandGap']
#                 color4FilledBandGap='b'
#                 if 'color4FilledBandGap' in kwargs: color4FilledBandGap=kwargs['color4FilledBandGap']
#             
#                 if isFilledBandGap:
#                     
#                     x=distances.tolist()+list(reversed(distances.tolist()))
#                     y=bands[:, upside[1]].tolist()+list(reversed(bands[:, downside[1]].tolist()))
#                 
#                     plt.fill(x, y, c=color4FilledBandGap)
#             elif btype[self.calculation] == 1:
#                 # up
#                 upside_up, downside_up, DirectBandgap_up=self.showBandInfo(bands[:,:,0], Efermi, tolerance4bandside=tolerance4bandside)
#                 # down
#                 upside_dn, downside_dn, DirectBandgap_dn=self.showBandInfo(bands[:,:,1], Efermi, tolerance4bandside=tolerance4bandside)
#                 bands -= np.min([downside_up[-1], downside_dn[-1]])
#                 for i in range(0, bands.shape[1]): plt.plot(distances, bands[:, i, 1], 'b', lw=1) # minority
#                 for i in range(0, bands.shape[1]): plt.plot(distances, bands[:, i, 0], 'r', lw=2) # majority
#         
#                 # fill bandagap
#                 # minority
#                 isFilledBandGap_dn=False
#                 if 'isFilledBandGap_dn' in kwargs: isFilledBandGap_dn=kwargs['isFilledBandGap_dn']
#                 color4FilledBandGap_dn='c'
#                 if 'color4FilledBandGap_dn' in kwargs: color4FilledBandGap_dn=kwargs['color4FilledBandGap_dn']
#             
#                 if isFilledBandGap_dn:
#                     x=distances.tolist()+list(reversed(distances.tolist()))
#                     y=bands[:, upside_dn[1], 0].tolist()+list(reversed(bands[:, downside_dn[1], 0].tolist()))
#                 
#                     plt.fill(x, y, c=color4FilledBandGap_dn)
#                 
#                 # majority
#                 isFilledBandGap_up=False
#                 if 'isFilledBandGap_up' in kwargs: isFilledBandGap_up=kwargs['isFilledBandGap_up']
#                 color4FilledBandGap_up='m'
#                 if 'color4FilledBandGap_up' in kwargs: color4FilledBandGap_up=kwargs['color4FilledBandGap_up']
#                 
#                 if isFilledBandGap_up:
#                     x=distances.tolist()+list(reversed(distances.tolist()))
#                     y=bands[:, upside_up[1], 0].tolist()+list(reversed(bands[:, downside_up[1], 0].tolist()))
#                 
#                     plt.fill(x, y, c=color4FilledBandGap_up)
#                 
#         # plot high symmetry lines
#         hspointsDistances, hspointsSymbols=self._plotHighSymmetryLines(path_of_structrues, plt)
#         plt.axhline(y=0, linestyle="--", color='k') # fermi line
# 
#         plt.xlim(distances_of_all[0], distances_of_all[-1])
#         plt.ylim(erange[0],erange[1])
#         
#         # ylim=plt.ylim()
#         # deltay=ylim[1]-ylim[0]
#         
#         plt.tick_params(labelsize=18)
#         plt.xticks(hspointsDistances, hspointsSymbols, fontsize=18)
#         #plt.yticks(fontsize=18)
#         plt.ylabel('$\mathregular{E\ (ev)}$', fontsize=18).set_fontweight('bold')
#         
#         plt.tight_layout()
#         plt.savefig(filename+'.png', dpi=600)
#         plt.clf() # clear the figure
# =============================================================================
       
    def plotBand(self, path_of_structrues, erange, subplot=None, **kwargs):
        """
        plot band structure from EIGENVAL file.
        
        Arguments:
            path_of_structures:
            erange: range of plotted energy. e.g. [-4, 4]
            
            kwargs:
                filename: output filename. e.g. bandfat-Li
                isFilledBandGap (default=False): whether to fill the band gap.
                
                for non-spin:
                    isFilledBandGap (default: False):
                    color4FilledBandGap (default: 'b'): the color of filling band gap.
                for spin:
                    isFilledBandGap_up (default: False):
                    isFilledBandGap_dn (default: False):
                    color4FilledBandGap_up (default='m'): color of filling band gap.
                    color4FilledBandGap_dn (default='c'): color of filling band gap.
                    
                tolerance4bandside (default=1e-3):
                tolerance4kpoint (defult=1e-6):
        """
        
        
        btype=self.btype
        
        isSubplot=False if subplot is None else True
        if not(isSubplot):
            filename='band' if not('filename' in kwargs) else kwargs['filename']
        
        tolerance4bandside=1e-6
        if 'tolerance4bandside' in kwargs: tolerance4bandside=kwargs['tolerance4bandside']
# =============================================================================
#         tolerance4kpoint=1e-6
#         if 'tolerance4kpoint' in kwargs: tolerance4kpoint=kwargs['tolerance4kpoint']
# =============================================================================
        
        if not(isSubplot):
            fig=plt.figure(figsize=(8,8)) # size of canvas: 8x8 inch
        
            outer_grid = gridspec.GridSpec(1, 1)
            inner_grid = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer_grid[:, :], 
                                                          wspace=0.0, hspace=0.0)
        
            # band
            subplot = plt.Subplot(fig, inner_grid[0:1, 0:1])
        
        paths=[]
        if isinstance(path_of_structrues, str): 
            paths=[path_of_structrues]
        elif isinstance(path_of_structrues, list) or isinstance(path_of_structrues, np.ndarray):
            paths=path_of_structrues
        
        # prepare data
        distances_of_all=[]
        for ip in range(0, len(paths)):
            path=paths[ip]
            data=self.getReducedBand(path=path)
            kpoints=data['kpoints']
            bands=data['bands']
        
            distances=self.calcKpointsDistance(kpoints, path)
            
            # process distances
            if ip == 0:    
                distances_of_all=distances.tolist()
            else:
                distances += distances_of_all[-1]
                distances_of_all=distances_of_all+distances.tolist()
                    
            # plot band
            Efermi=self.getEfermi(path)
            if btype[self.calculation] == 0:
                upside, downside, DirectBandgap=self.showBandInfo(bands, Efermi, tolerance4bandside=tolerance4bandside)
                # move Fermi to VBM
                bands -= downside[-1]
                for i in range(0, bands.shape[1]): subplot.plot(distances, bands[:, i], 'k', lw=1)
            
                # fill bandagap
                isFilledBandGap=False
                if 'isFilledBandGap' in kwargs: isFilledBandGap=kwargs['isFilledBandGap']
                color4FilledBandGap='b'
                if 'color4FilledBandGap' in kwargs: color4FilledBandGap=kwargs['color4FilledBandGap']
            
                if isFilledBandGap:
                    
                    x=distances.tolist()+list(reversed(distances.tolist()))
                    y=bands[:, upside[1]].tolist()+list(reversed(bands[:, downside[1]].tolist()))
                
                    subplot.fill(x, y, c=color4FilledBandGap)
            elif btype[self.calculation] == 1:
                # up
                upside_up, downside_up, DirectBandgap_up=self.showBandInfo(bands[:,:,0], Efermi, tolerance4bandside=tolerance4bandside)
                # down
                upside_dn, downside_dn, DirectBandgap_dn=self.showBandInfo(bands[:,:,1], Efermi, tolerance4bandside=tolerance4bandside)
                bands -= np.min([downside_up[-1], downside_dn[-1]])
                for i in range(0, bands.shape[1]): subplot.plot(distances, bands[:, i, 1], 'b', lw=1) # minority
                for i in range(0, bands.shape[1]): subplot.plot(distances, bands[:, i, 0], 'r', lw=2) # majority
        
                # fill bandagap
                # minority
                isFilledBandGap_dn=False
                if 'isFilledBandGap_dn' in kwargs: isFilledBandGap_dn=kwargs['isFilledBandGap_dn']
                color4FilledBandGap_dn='c'
                if 'color4FilledBandGap_dn' in kwargs: color4FilledBandGap_dn=kwargs['color4FilledBandGap_dn']
            
                if isFilledBandGap_dn:
                    x=distances.tolist()+list(reversed(distances.tolist()))
                    y=bands[:, upside_dn[1], 0].tolist()+list(reversed(bands[:, downside_dn[1], 0].tolist()))
                
                    subplot.fill(x, y, c=color4FilledBandGap_dn)
                
                # majority
                isFilledBandGap_up=False
                if 'isFilledBandGap_up' in kwargs: isFilledBandGap_up=kwargs['isFilledBandGap_up']
                color4FilledBandGap_up='m'
                if 'color4FilledBandGap_up' in kwargs: color4FilledBandGap_up=kwargs['color4FilledBandGap_up']
                
                if isFilledBandGap_up:
                    x=distances.tolist()+list(reversed(distances.tolist()))
                    y=bands[:, upside_up[1], 0].tolist()+list(reversed(bands[:, downside_up[1], 0].tolist()))
                
                    subplot.fill(x, y, c=color4FilledBandGap_up)
                
        # plot high symmetry lines
        hspointsDistances, hspointsSymbols=self._plotHighSymmetryLines(path_of_structrues, subplot)
        subplot.axhline(y=0, linestyle="--", color='k') # fermi line

        subplot.set_xlim(distances_of_all[0], distances_of_all[-1])
        subplot.set_ylim(erange[0],erange[1])
        
        # ylim=subplot.ylim()
        # deltay=ylim[1]-ylim[0]
        
        subplot.tick_params(labelsize=18)
        #subplot.xaxis.set_ticks(hspointsDistances, hspointsSymbols)#, fontsize=18)
        #subplot.set_xticks(hspointsDistances, list(hspointsSymbols))#, fontsize=18)
        #plt.xticks(hspointsDistances, hspointsSymbols, fontsize=18)
        locs = subplot.set_xticks(hspointsDistances, minor=False)
        labels = subplot.set_xticklabels(hspointsSymbols, minor=False)
        #subplot.yticks(fontsize=18)
        subplot.set_ylabel('$\mathregular{E\ (ev)}$', fontsize=18).set_fontweight('bold')
        
        if not(isSubplot):
            fig.add_subplot(subplot)
            
            fig.tight_layout()
            fig.savefig(filename+'.png', dpi=600)
            fig.clf() # clear the figure
        
    def plotScatterProjectedBand(self, projectedElements, filename, erange, scale, path, **kwargs):
        """
        plot scatter-style projected band. Note that it'll plot the scatters of given elements in the order from big to small on each (nb, nk).
        
        Arugments:
            projectedElements: python list of two projected elements. [(element, [norbit, color]),...] e.g. [('Li', [1, 'r']), ('P', [1,'g'])]
            filename: output filename. e.g. bandfat-Li
            erange: range of plotted energy. e.g. [-4, 4]
            scale: scale of projected scatter (circle).
            
            kwargs:
                For SOC calculation, there is a additional parameter:
                    m (default=0): the range of value is [0, 4]. 0: m_toatal; 1: m_x; 2: m_y; 3: m_z.
                
                tolerance4bandside (default=1e-3)
        """
        pass
        
    def plotCustomScatterProjectedBand(self, projectedElements, filename, erange, scales, path, **kwargs):
        """
        plot custom scatter-style projected band. Note that it'll plot the scatters of given elements in the order from big to small on each (nb, nk).
        
        Arugments:
            projectedElements: python list of two projected elements. [(element, [norbit, color]),...] e.g. [('Li', [1, 'r']), ('P', [1,'g'])]
            filename: output filename. e.g. bandfat-Li
            erange: range of plotted energy. e.g. [-4, 4]
            scales: custom scales of projected scatter (circle). The valid format:
                [[nb, nk0, nk1, scale], [nb, nk0, nk1, scale],..., [-1, -1, -1, scale]]
                         case1                   case2                     other
            kwargs:
                For SOC calculation, there is a additional parameter:
                    m (default=0): the range of value is [0, 4]. 0: m_toatal; 1: m_x; 2: m_y; 3: m_z.
                    
                tolerance4bandside (default=1e-3)
        """
        def getScale(nb, nk, scales):
            """
            get the scale value for given (nb, nk)
            """
            scales=np.array(scales)
            
            scale=None
            for i in range(0, scales.shape[0]-1):    
                if nb == scales[i][0] and (nk >= scales[i][1] and nk <= scales[i][2]):
                    scale=scales[i][3]
            if scale is None:
                scale=scales[-1][3]
                
            return scale
        
        pass

    def plotGradientBand(self, projectedElements, filename, erange, path, **kwargs):
        """
        plot the gradient band.
        
        Arguments:
            projectedElements: python list of two projected elements. [(element, [norbit, color]),...] e.g. [('Li', [1, 'r']), ('P', [1,'g'])]
            filename: output filename. e.g. bandfat-Li
            erange: range of plotted energy. e.g. [-4, 4]
            
            kwargs:
                For SOC calculation, there is a additional parameter:
                    m (default=0): the range of value is [0, 4]. 0: m_toatal; 1: m_x; 2: m_y; 3: m_z.
                    
                tolerance4bandside (default=1e-3)
        """
        pass
    
    def plotPDOS(self, path, erange, subplot=None, **kwargs):
        """
        Args:
            path:
            erange:
        
            kwargs:
                filename:
                isVertical (default=False)
                hasTotal (default=True)
        """
        def getMaxValue(totalDOS, energies, erange):
            maxv=None
            for i in range(0, len(energies)):
                if (energies[i] >= erange[0]) and (energies[i] < erange[1]):
                    if (maxv is None) or (totalDOS[i] > maxv): maxv=totalDOS[i]
            return maxv
        
        btype=self.btype
        
        isSubplot=False if subplot is None else True
        if not(isSubplot):
            filename='pdos' if not('filename' in kwargs) else kwargs['filename']
        isVertical=False if not('isVertical' in kwargs) else kwargs['isVertical']
        hasTotal=True if not('hasTotal' in kwargs) else kwargs['hasTotal']
        
        if not(isSubplot):
            fig=plt.figure(figsize=(8,8)) # size of canvas: 8x8 inch
            outer_grid = gridspec.GridSpec(1, 1)
            inner_grid = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer_grid[:, :], 
                                                      wspace=0.0, hspace=0.0)
        
            # pdos
            subplot = plt.Subplot(fig, inner_grid[0:1, 0:1])

        a, b, volume=self.getLatticeVector(path=path)
        
        data=self.getDOS(path=path)
        energies=data['energies']
        totalDOS=data['totalDOS']
        cumulativeDOS=data['cumulativeDOS']
        atoms=data['atoms']
        elements=data['elements']
        
        
        for element, pdos in elements.items():
            if not(isVertical):
                if btype[self.calculation] == 0:
                    subplot.plot(energies, pdos.sum(axis=1)/volume, lw=1, label=element)
                elif btype[self.calculation] == 1:
                    subplot.plot(energies, pdos[:,0:-1:2].sum(axis=1)/volume, lw=1, label=element)
                    subplot.plot(energies, -pdos[:,1:-1:2].sum(axis=1)/volume, lw=1)
            else:
                if btype[self.calculation] == 0:
                    subplot.plot(pdos.sum(axis=1)/volume, energies, lw=1, label=element)
                elif btype[self.calculation] == 1:
                    subplot.plot(pdos[:,0:-1:2].sum(axis=1)/volume, energies, lw=1, label=element)
                    subplot.plot(-pdos[:,1:-1:2].sum(axis=1)/volume, energies, lw=1)
        
        if hasTotal:
            if not(isVertical):
                if btype[self.calculation] == 0:
                    subplot.plot(energies, totalDOS/volume, 'k', lw=1, label='total')
                elif btype[self.calculation] == 1:
                    subplot.plot(energies, totalDOS[:, 0]/volume, 'k', lw=1, label='total')
                    subplot.plot(energies, -totalDOS[:, 1]/volume, 'k', lw=1)
            else:
                if btype[self.calculation] == 0:
                    subplot.plot(totalDOS/volume, energies, 'k', lw=1, label='total')
                elif btype[self.calculation] == 1:
                    subplot.plot(totalDOS[:, 0]/volume, energies, 'k', lw=1, label='total')
                    subplot.plot(-totalDOS[:, 1]/volume, energies, 'k', lw=1)
        
        # modify x and y ticks
        if not(isVertical):
            subplot.set_xlim(erange[0],erange[1])
            # yticks
            if btype[self.calculation] == 0:
                subplot.set_ylim(0, getMaxValue(totalDOS=totalDOS, energies=energies, erange=erange)/volume)
            elif btype[self.calculation] == 1:
                maxv=np.max([getMaxValue(totalDOS=totalDOS[:, 0], energies=energies, erange=erange), 
                             getMaxValue(totalDOS=totalDOS[:, 1], energies=energies, erange=erange)])
                subplot.set_ylim(-maxv/volume*1.1, maxv/volume*1.1)
            subplot.set_xlabel('$\mathregular{E\ (ev)}$', fontsize=18).set_fontweight('bold')
            subplot.set_ylabel('$\mathregular{PDOS\ (\AA^{-3})}$', fontsize=18).set_fontweight('bold')
        else:
            subplot.set_ylim(erange[0],erange[1])
            # xticks
            if btype[self.calculation] == 0:
                subplot.set_xlim(0, getMaxValue(totalDOS=totalDOS, energies=energies, erange=erange)/volume)
            elif btype[self.calculation] == 1:
                maxv=np.max([getMaxValue(totalDOS=totalDOS[:, 0], energies=energies, erange=erange), 
                             getMaxValue(totalDOS=totalDOS[:, 1], energies=energies, erange=erange)])
                subplot.set_xlim(-maxv/volume*1.1, maxv/volume*1.1)
            subplot.set_ylabel('$\mathregular{E\ (ev)}$', fontsize=18).set_fontweight('bold')
            subplot.set_xlabel('$\mathregular{PDOS\ (\AA^{-3})}$', fontsize=18).set_fontweight('bold')
        
        subplot.tick_params(labelsize=18)
        # subplot.xticks(hspointsDistances, hspointsSymbols, fontsize=18)
        #subplot.yticks(fontsize=18)
        # subplot.xlabel('$\mathregular{E\ (ev)}$', fontsize=18).set_fontweight('bold')
        # subplot.ylabel('$\mathregular{PDOS\ (\AA^{-3})}$', fontsize=18).set_fontweight('bold')
        
        subplot.legend(loc=1, numpoints=2,
                   prop=(FontProperties(weight="bold", size=18)), frameon=True)
        
        if not(isSubplot):
            fig.add_subplot(subplot)
            
            fig.tight_layout()
            fig.savefig(filename+'.png', dpi=600)
            fig.clf() # clear the figure
        

    def plotBand_PDOS(self, path_of_bands, path_of_dos, erange, filename='band2pdos', **kwargs):
        """
        """
        fig=plt.figure(figsize=(8,8)) # size of canvas: 8x8 inch
        
        outer_grid = gridspec.GridSpec(1, 1)
        inner_grid = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_grid[:, :], 
                                                          wspace=0.05, hspace=0.0)
        
        # band
        subplot0 = plt.Subplot(fig, inner_grid[0:1, 0:2])
        
        self.plotBand(path_of_structrues=path_of_bands, erange=erange, subplot=subplot0)
        
        fig.add_subplot(subplot0)
        
        # band
        subplot2 = plt.Subplot(fig, inner_grid[0:1, 2:3])
        
        self.plotPDOS(path=path_of_dos, erange=erange, subplot=subplot2, isVertical=True, hasTotal=True)
        subplot2.set_ylabel('')
        subplot2.set_yticks([])
        subplot2.set_xticks([])
        
        fig.add_subplot(subplot2)
        
        fig.tight_layout()
        fig.savefig(filename+'.png', dpi=600)
        fig.clf() # clear the figure
   
    
   
    
# -------------------- test --------------------

# b.plotBand(path_of_structrues=[os.getcwd(), os.getcwd()], erange=[-4, 10], filename='band', isFilledBandGap_up=True)
# b.plotPDOS(path=os.getcwd(), erange=[-10, 10], filename='pdos')

# b.plotBand(path_of_structrues=[os.getcwd(), os.getcwd()], erange=[-4, 10], filename='band', isFilledBandGap_up=True)
# b.plotPDOS(path=os.getcwd(), erange=[-10, 10], filename='pdos', isVertical=True)

b=PlotBand()

if b.Band:
    print('band: {}'.format(b.path_of_band))
    print('erange: {}'.format(b.erange))
    
    b.plotBand(path_of_structrues=b.path_of_band,  
               erange=b.erange)

if b.PDOS:
    print('dos: {}'.format(b.path_of_dos))
    print('erange: {}'.format(b.erange))
    
    b.plotPDOS(path=b.path_of_dos, 
               erange=b.erange)

if b.Band_PDOS:
    print('band: {}'.format(b.path_of_band))
    print('dos: {}'.format(b.path_of_dos))
    print('erange: {}'.format(b.erange))
    
    b.plotBand_PDOS(path_of_bands=b.path_of_band, 
                    path_of_dos=b.path_of_dos, 
                    erange=b.erange)
