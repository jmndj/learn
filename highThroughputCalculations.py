#!/usr/bin/env python
'''
Created on Nov 3, 2022

@author: fu
'''
import os
from time import sleep

class HighThroughputCalculations():
    def __init__(self, maxJobs, type_of_system, path=None):
        """Arguments:
            maxJobs:
            type_of_system: 'pbs' or 'lsf'
            threshold_of_CE:
        
        """
        self.maxJobs=maxJobs
        self.type_of_system=type_of_system
        
        if path is  None:
            self.path=os.getcwd()
        else:
            self.path=path
    
    def replaceElements(self,composition, path_of_orig, path_of_output):
        """
        composition: {'Na': ['Na', 'K', 'Rb', 'Cs'],
                      'Cl': ['Cl', 'Br', 'I']}
        path_of_orig:
        path_of_output:
        """
        os.mkdir(path_of_output)
        m = 0
        #pass
        with open('%s/0-NaCl.vasp'%(path_of_orig), 'r') as file:
            lines = file.readlines()
        for i in composition['Na']:
            for j in composition['Cl']:
                n = 5
                lines[n] = lines[n].replace("Cl", "{}".format(j)) 
                lines[n] = lines[n].replace("Na", "{}".format(i))
                with open('%s/%s-%s%s.vasp'%(path_of_output,m,i,j), 'w+') as file:
                    file.writelines(lines)
                m = m+1

        
    
            
    def generateCalculatedDirectory(self, path):
        """
        """
        structures=sorted([s for s in os.listdir(path) if s.endswith('.vasp')], key=lambda s: int(s.split('-')[0]))
        for s in structures:
            destinaiton_path=path+'/%s' %s.split('-')[0]
            if not(os.path.exists(destinaiton_path)):
                os.mkdir(path+'/%s' %s.split('-')[0])
                os.system('cp %s/%s %s/POSCAR' %(path, s, destinaiton_path))
                
                os.system('touch %s/wait' %destinaiton_path)
                
    def getList(self, type, path):
        """
        Arguments:
            type: 'wait' or 'calculate' or 'error'
        """
        print('a', os.getcwd(), path)
        structures=sorted([s for s in os.listdir(path) if os.path.isdir(path+'/%s' %s)], key=lambda s: int(s))
        
        slist=[]
        for s in structures:
            if os.path.exists(path+'/%s/%s' %(s, type)):
                slist.append(s)
        
        return slist
    
    def check(self, filename='run.log'):
        """
        check the given calculation
        """
        string=os.popen('grep "reached required accuracy" %s' %filename).readline()
        if string != '':
            return True
        else:
            return False
        
    def isContinue(self, path):
        """
        """
        isContinue=True
        if os.path.exists(path+'/over'):
            isContinue=False

        return isContinue
    
    def prepareSubmit(self, id, pref, filename='submit.sh'):
        """
        modify the name of jobs: 'pref-id', e.g., 'yhfu-1'.
        
        Arguments:
            id:
            pref:
            filename:
        """
        os.system('cp ../%s .' %(filename))

        with open(filename, 'r') as file:
            lines = file.readlines()

        # 修改列表中的第n行（注意Python的索引是从0开始的）
        n = 1
        lines[n] = lines[n].replace("XXX", "yhfu-{}".format(id))      
        with open(filename, 'w') as file:
            file.writelines(lines)
        
    def log(self, string, path, filename='monitor.log'):
        """
        """
        outfile=open('%s/%s' %(path, filename), 'a')
        outfile.write('%s\n' %string)
        outfile.close()
        
    def run(self, pref_jobs, path):
        """
        """
        
        def check_finished():
            # check finished
            finished=self.getList(type='finished', path=path)
            error=self.getList(type='error', path=path)
            print('%s  finished: %d' %(filename, len(finished)), finished)
            print('%s  error: %d' %(filename, len(error)), error)
            self.log('%s  finished: %d %s\n' %(filename, len(finished), '['+', '.join(finished)+']'), self.path)
            self.log('%s  error: %d %s\n\n' %(filename, len(error), '['+', '.join(error)+']'), self.path)
            for f in finished:
                print('f->', f)
                os.chdir('%s/%s' %(path, f))
                if self.check():
                    os.system("grep 'energy  without entropy' OUTCAR | tail -1 | awk '{print $4}' > energy")
                    os.system('rm calculate wait')
                else:
                    os.system('rm finished calculate wait; touch error')
                os.chdir(path)
                
        self.generateCalculatedDirectory(path)
        os.chdir(path)

        filename=os.path.basename(path)
    
        # monitor
        while self.isContinue(path=self.path):
            # keep jobs to maximum
            waits=self.getList(type='wait', path=path)

            # update jobs on cluster
            calculates=self.getList(type='calculate', path=path)
            waits=self.getList(type='wait', path=path)
            prepares=sorted(list(set(waits)-set(calculates)), key=lambda s:int(s))
            print('%s  calculates: %d' %(filename, len(calculates)), calculates)
            print('%s  waits: %d' %(filename, len(waits)), waits)
            print('%s  prepares: %d' %(filename, len(prepares)), prepares)
            self.log('%s  calculates: %d %s\n' %(filename, len(calculates), '['+', '.join(calculates)+']'), self.path)
            self.log('%s  waits: %d %s\n' %(filename, len(waits), '['+', '.join(waits)+']'), self.path)
            self.log('%s  prepares: %d %s\n' %(filename, len(prepares), '['+', '.join(prepares)+']'), self.path)
          
            if prepares == []:
                check_finished()
                break
            
            if len(calculates) < self.maxJobs:
                p=prepares[0]
                print('p->', p)
                os.chdir('%s/%s' %(path, p))
                print(os.getcwd())
                self.prepareSubmit(id=p, pref=pref_jobs)
                if self.type_of_system.lower() == 'pbs':
                    os.system('cp ../INCAR_* .; qsub submit.sh')
                elif self.type_of_system.lower() == 'lsf':
                    os.system('cp ../INCAR_* .; bsub < submit.sh')
                os.system('touch calculate; rm wait')
                os.chdir(path)
                
            check_finished()
            
            sleep(30)
        os.system('rm over')
        
    def batch(self, orig_path, filename='dp'):
        from itertools import combinations
        
        
        if not(self.isContinue(path=self.path)):
            exit()

        # filename='dp'
            
        path=self.path+'/'+filename
        print(path)
        if not(os.path.isdir(path)):
            os.mkdir(path)
        os.chdir(path)
        print(os.getcwd())
            
        os.system('cp %s/INCAR_* %s/submit.sh %s/run.sh  %s' %(self.path, self.path, self.path, path))
        os.system('cp %s/*.vasp %s' %(orig_path, path))
        # os.system('cp %s/OPTCELL-* %s' %(orig_path, path))
        
        self.run(pref_jobs='%s' %('dp'), path=path)
                
        os.chdir(self.path)
        os.system('rm over')
        
        
# ---------- test ----------
path=os.getcwd()
r=HighThroughputCalculations(maxJobs=2, type_of_system='pbs', path=path)
composition={'Na': ['Na', 'K', 'Rb', 'Cs'],
              'Cl': ['Cl', 'Br', 'I']}
r.replaceElements(composition,path_of_orig = path+'/orig/calculate',path_of_output='orig/calculate_opt')
r.batch(orig_path=path+'/orig/calculate_opt')
