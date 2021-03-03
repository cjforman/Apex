'''
Created on 26 Feb 2021

@author: Chris
'''
import sys
from FileIO import FileIO

class PDBIO(FileIO):
    '''
    A FileIO object that handles file operations on PDBs 
    '''

    def __init__(self, params):
        '''
        Constructor
        '''
        self.params = params

    # takes an input pdb file and strips away the ".pdb" header leaving just the fileroot
    def fileRootFromFilename(self, pdbFilename):
        return pdbFilename.replace('.pdb', '', 1)
        
    def lineFromPDBAtom(self, atom):
        try:
            l='ATOM {: >06d} {: <4}{:1}{:3} {:1}{: >4d}{:1}   {: >8.3f}{: >8.3f}{: >8.3f}{: >6.2f}{: >6.2f}      {: <4}{: >2}{: >2}\n'.format(int(atom[0]), atom[1], atom[2], atom[3], atom[4], int(atom[5]), atom[6], float(atom[7]), float(atom[8]), float(atom[9]), float(atom[10]), float(atom[11]), atom[12],atom[13],atom[14])
        except:
            print("unable to write atom to string: ")
            print(atom)
            sys.exit()
        return l
 
    def readPDB(self, filename, atomTypes, heartBeat=0):
        return self.parsePDB(self.loadData(filename, 'txt'), atomTypes, heartBeat=heartBeat)
    
    def parsePDB(self, lines, atomTypes=['ATOM'], heartBeat=0):
        #parse data in lines  
        atoms = []
        for count, line in enumerate(lines):
            if not heartBeat==0 and count % heartBeat == 0:
                print("Parsing Atom: ", count) 
            
            parse = False
            for entry in atomTypes:
                if line[0:4] == entry:
                    parse = True
                    break
            if parse:
                try:
                    a = self.parsePDBLine(line)
                    if not a in atoms:
                        atoms.append(a)
                except:
                    print("line: " + line + " not understood")
                    sys.exit()
        return atoms
    
    
    def parsePDBLine(self, line):
        l=line
        atom_seri = int(l[6:11])
        atom_name = l[12:16].split()[0]
        alte_loca = l[16]
        resi_name = l[17:20].split()[0]
        chai_iden = l[21]
        resi_numb = int(l[22:26])
        code_inse = l[26]
        atom_xcoo = float(l[30:38])
        atom_ycoo = float(l[38:46])
        atom_zcoo = float(l[46:54])
        try:
            atom_occu = float(l[54:60])
        except:
            atom_occu=0.0
    
        try:
            atom_bfac = float(l[60:66])
        except:
            atom_bfac=0.0    
        
        try:
            seg_id = l[72:76]
        except:
            seg_id=' '
    
        try:
            atom_symb = l[76:78].split()[0]
        except:
            try:
                atom_symb = l[68]
            except:
                atom_symb= ' '
    
        try:
            charge=l[78:80]
        except:
            charge=' '
    
        return [atom_seri, atom_name, alte_loca, resi_name, chai_iden, resi_numb, code_inse, atom_xcoo, atom_ycoo, atom_zcoo, atom_occu, atom_bfac,seg_id,atom_symb,charge]