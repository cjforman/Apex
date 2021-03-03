'''
Created on 2 March 2021

@author: Chris
'''

from project import Project
from PDBProc import PDBProc as PDBProc
import sys

class PDBAnalyzer(Project):

    def __init__(self, settings):
        # copy any settings over
        self.settings = settings
        
    def doPDBProc(self):
        PDBProcObj = PDBProc(self.settings)
        PDBProcObj.doAnalysePDB()
                
         
if __name__ == '__main__':
    
    config_filename = sys.argv[1]
    
    PDBAnalyzer_obj = PDBAnalyzer.from_configfile(config_filename)
    
    # perform tasks using analyser
    PDBAnalyzer_obj.doPDBProc()
    
    print("Done")         
    
