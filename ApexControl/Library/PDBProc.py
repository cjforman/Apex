'''
Created on 26 Feb 2021

@author: Chris
'''
from PDBIO import PDBIO
from PlotData import plotData as plotData
import sys
import numpy as np


class PDBProc(PDBIO):
    '''
    A FileIO object that handles file operations on PDBs 
    '''

    def __init__(self, params):
        self.usage ='''
        The config file controls the PDBProc Code. 
        
        Mandatory:
        
        "PDBFilename": <pdbfilename>  in format: "fileroot.pdb" 
        
        "command": <command>
        
        "atomTypes": ["ATOM", "HETATM"]   can choose to ignore HETATM by leaving out HETATM from this list. 

        OptionaL:
        
        'outputFileroot': "outfileroot"  defines the fileroot for the output. Different commands append their own header to create output filename.  
                                         If outputFileroot is not specified then the fileroot from the original pdb is used as the output filename.    
    
        'verbose': [0|1]  can be 0 or 1. If 1 prints out full results on screen. Omission defaults to 0
        
        'heartBeat': period  If >0 prints out results on screen every period entries to reassure user the 
                             process is still active on long operations. Omission defaults to 0.


        Example:
          {
            "settings": {
                    "PDBFilename": "MASP2.pdb",
                    "command": "resInfo",
                    "atomTypes": ["ATOM"],
                    "verbose": 0,
                    "heartBeat": 1000,
                    "outputFileroot": "MASP2",
                    "resInfoType": ["torsion"],
                    "plotData": ["Ramachandran"],
                    "plotDict": {"PRO":"g+", "TYR": "bx", "DEF": "rx"},
                    "Grid": 1,
                    "phiMin": -180, 
                    "phiMax":  180, 
                    "psiMin": -180, 
                    "psiMax":  180,
                    "Title": "MASP2 Ramachandran Plot"
                }
        }

        
        'command' can be one of:
        
        'resInfo':   Provides a structure of the PDB atom info as a json with filename 'fileroot.res'
                     Saves the Settings config file in that structure for plotting control with plotData. 
                     
                    Required parameters:
                    'resInfoType': ['torsion']  controls which residue information is including in the json out put file.
                    
                    e.g. 'torsion': returns the phi, psi and omega values for each residue in the chain. 
                                    Return None for the end Phi, Psi and omegas.
                                   
                                    Also creates a csv text file called: fileroot.tor

                                    and checks for a plotData style called "ramachandran" which creates a ramaChandran plot.

        
        'plotData':  ['ramachandran']
                     invokes a plotData command with each of the specified plotting routine keywords. 
                     Call each plot data routines listed and uses the json of the output data from the command as input data for that plot.
                     It's also possible to call that plotData from the command line using the saved json file as input.  
                     Settings for the plotting are included in the config json.

                    "plotDict": {"PRO":"g+", "TYR": "bx", "DEF": "rx"},  specifies the plotting symbol for each residue.
                    
                    The following commands are self=explanatory
                    "legend":
                    "Grid": 1,
                    "phiMin": -180, 
                    "phiMax":  180, 
                    "psiMin": -180, 
                    "psiMax":  180,
                    "Title": "MASP2 Ramachandran Plot"

        
            
        '''
        self.params = params
        
    # function makes sure that the config file makes sense for the given PDB processing command before calling that command
    def doAnalysePDB(self):
        
        # make sure there is PDBFilename
        try:
            self.PDBFilename = self.params['PDBFilename']
        except KeyError as e:
            print("Must define 'PDBFilename' in the config file to process a PDB.", e)
            print(self.usage)
            sys.exit()
        
        # make sure that there is a command to perform on the filename
        try:
            command = self.params['command']
        except KeyError as e:
            print("Must provide a command to perform on the PDB file.", e)
            print(self.usage)
            sys.exit()

        try:
            self.outputFileroot = self.params['outputFileroot']
        except KeyError:
            self.outputFileroot = self.fileRootFromFilename(self.PDBFilename)

        try:
            self.verbose = self.params['verbose']
        except KeyError:
            self.verbose = 0

        try:
            self.heartBeat = self.params['heartBeat']
        except KeyError:
            self.heartBeat = 0

                    
        print("Performing ", command, " on ", self.PDBFilename, " with outputfileroot: ", self.outputFileroot, " with heartBeat:", self.heartBeat)

        
        # res info command gathers info about the residues of the PDB such as torsion 
        if command=='resInfo':
            try:
                self.resInfoType = self.params['resInfoType']
            except KeyError as e:
                print("Command: resInfo requires a list of 'resInfoTypes', such as: 'torsion'. Key Error", e)
                print(self.usage)
                sys.exit()

            results = self.analyseResidues()
    
    
        # save the json of the return data file from the command  
        self.saveData(self.outputFileroot + ".json", results)

        # calls each plot function defined confiog file
        for plotCommand in self.params['plotData']:
            plotData(self.outputFileroot + ".json", plotCommand)

    
        return results
    
    # seperates a list of atoms into chains based on the chain letter.
    # if chains letters are not defined then all atoms will go into the same chain.
    # if the chains are not all the same length then the program will end
    def AtomsToChains(self, atoms):
        ''' Breaks a list of atoms into chains based on the chain letters. 
            If there are multiple stuctures the chains from each structure will all be 
            assigned to the same dictionary. May need to change that in the future. '''
        
        # set up a dictionary to be the output set of chains indexed by chain letter
        chains = {}
        # loop through the atoms and append each atom to the chain to which it belongs
        # creates a chain if it is the first time a chain letter is encountered
        for count, atom in enumerate(atoms):
            try:
                chains[atom[4]].append(atom)
            except KeyError:
                chains[atom[4]] = [atom]

            if not self.heartBeat==0 and count % self.heartBeat == 0:
                print("Parsing Atom: ", count, " into chain: ", atom[4]) 
                
        if self.verbose:
            print(str(len(chains)) + " chain(s) found: id(s): ", [chain for chain in chains])
        
        chainLengths = [len(chains[chain]) for chain in chains ]
        
        if self.verbose:
            print("Chain lengths are: ", chainLengths)
        
        # if the set of chain lengths does not have one and only number in it then the chain lengths are different lengths    
        if len(set(chainLengths))!=1:
            print("Error: chains are not the same length in file:", self.PDBFilename)
            sys.exit()
            
        return chains

    # creates a data structure which is a dictionary of dictionaries indexed by chain letter and resIDNumbers and puts all the atoms in each residues in each subdictionary for that residue.
    # also creates a list of resIDs as strings and stores it in each chain dictionary 
    def generateResInfoStructureFromChains(self, chains):
        # create output dictionary
        resInfo = {}

        # go through the chain appending each atom to dictionary keyed under it's residue number and chain id.
        for chain in chains:
            # create sub dictionary for the chain
            resInfo[chain] = {}
            
            # add a list to contain the residue IDS in this chain.
            resInfo[chain]['resList'] = []
            
            # loop through each atom in the chain and add it to the appropriate residue ID. 
            for count, atom in enumerate(chains[chain]):
                try:
                    # if residue ID is already there and append this atom to the list for that residue
                    resInfo[chain][atom[5]]['atomList'].append(atom)
                    if not self.heartBeat==0 and count % self.heartBeat == 0:
                        print("Parsing Atom: ", count, " into chain: ", atom[4], " and residue: ", atom[5]) 
 
                except KeyError:
                    # if resid is not there then create the residue dictionary, 
                    # and start the atomList for that residue off.
                    # Also add the residue number to the list of residue keys in the chain in 
                    # the order they appear in the PDB 
                    resInfo[chain][atom[5]] ={}
                    resInfo[chain][atom[5]]['atomList'] = [atom]
                    resInfo[chain]['resList'].append(atom[5])
                    
        return resInfo
        
        
    def analyseResidues(self):
        ''' Analyses each residue and compiles lists based on information requested in the params'''
        
        try:
            atoms = self.readPDB(self.PDBFilename, atomTypes = self.params['atomTypes'], heartBeat=self.heartBeat)
        except KeyError as e:
            print("KeyError: ", e)
            sys.exit()
        
        # create a dictionary with all the chains
        chains = self.AtomsToChains(atoms)
        
        # generate the residue info structure
        resInfo = self.generateResInfoStructureFromChains(chains)

        # check to see which information is requested in the config file. 
        if 'torsion' in self.resInfoType:              
            # output list
            torsionList = [] 
            
            # if torsion is requested then loop through the chains
            for chain in chains:
                
                # extract the reslist foir the current chain to simplify code a little. 
                resList = resInfo[chain]['resList']

                # Deal with the first residues in the chain 
                resInfo[chain][resList[0]]['torsion'] = self.computeTorsions(None, 
                                                                             resInfo[chain][resList[0]]['atomList'], 
                                                                             resInfo[chain][resList[1]]['atomList'])
                
                # torsion info for first residue to output list
                torsionList.append([chain, 
                                    resList[0], 
                                    resInfo[chain][resList[0]]['torsion'][0], 
                                    resInfo[chain][resList[0]]['torsion'][1], 
                                    resInfo[chain][resList[0]]['torsion'][2]])
                count= 0 
                # compute torsions for all the middle residues. 
                for prevR, curR, nextR in zip(resInfo[chain]['resList'][0:-2], 
                                              resInfo[chain]['resList'][1:-1], 
                                              resInfo[chain]['resList'][2:]): 
                    
                    if not self.heartBeat==0 and count % self.heartBeat==0:
                        print("Computing torsion for chain: ", chain, " and residue: ", curR) 
                    
                    # compute the psi, phi and omega torsion angles for each residue - needs info from prev and next residues.
                    # assumes the order that residues come in the in List defines the neighbouring residues in the structure. 
                    resInfo[chain][curR]['torsion'] = self.computeTorsions(resInfo[chain][prevR]['atomList'], 
                                                                           resInfo[chain][curR]['atomList'], 
                                                                           resInfo[chain][nextR]['atomList'])
                
                    # append torsion info to the output list
                    torsionList.append([chain, 
                                        curR, 
                                        resInfo[chain][curR]['torsion'][0], 
                                        resInfo[chain][curR]['torsion'][1], 
                                        resInfo[chain][curR]['torsion'][2]])

                    count += 1

                # deal with the last residue in the chain
                resInfo[chain][resList[-1]]['torsion'] = self.computeTorsions(resInfo[chain][resList[-2]]['atomList'], 
                                                                              resInfo[chain][resList[-1]]['atomList'], 
                                                                              None)
                # output the last residue torsion info to the list
                torsionList.append([chain, 
                               resList[-1], 
                               resInfo[chain][resList[-1]]['torsion'][0], 
                               resInfo[chain][resList[-1]]['torsion'][1], 
                               resInfo[chain][resList[-1]]['torsion'][2]])


            if self.params['verbose']==1:
                [ print(tor) for tor in torsionList]

            # output the list of data for all chains to a csv file with a .tor header
            self.saveData(self.outputFileroot + '.tor', torsionList, delimiter=',')

        # create an output data structure to save which includes the data and the settings
        outData = {}
        outData['settings'] = self.params
        outData['resInfo'] = resInfo

        return outData
    
    def computeTorsions(self, res0, res1, res2):
        # find the necessary atoms in the residues to compute the desired dihedrals
        # if any are not found then set the relevant dihedral to None.
        try:
            #extract the relevant information from the residues and present as position vectors of atoms.
            C0 = [np.array([atom[7],atom[8],atom[9]]) for atom in res0 if atom[1]=='C'][0]
        except:
            phi = [None, 0]

        try:
            N1 = [np.array([atom[7],atom[8],atom[9]]) for atom in res1 if atom[1]=='N'][0]
        except:
            phi = [None, 0]
            psi = [None, 0]
        
        try:
            CA1 = [np.array([atom[7],atom[8],atom[9]]) for atom in res1 if atom[1]=='CA'][0]
        except:
            phi = [None, 0]
            psi = [None, 0]
            omega = [None, 0]

        try:
            C1 = [np.array([atom[7],atom[8],atom[9]]) for atom in res1 if atom[1]=='C'][0]            
        except:
            phi = [None, 0]
            psi = [None, 0]
            omega = [None, 0]
        
        try:
            N2 = [np.array([atom[7],atom[8],atom[9]]) for atom in res2 if atom[1]=='N'][0]
            
        except:
            psi = [None, 0]
            omega = [None, 0]

        try:
            CA2 = [np.array([atom[7],atom[8],atom[9]]) for atom in res2 if atom[1]=='CA'][0]
        except:
            omega = [None, 0]
            
        
        #compute phi; the C0'-N1-Ca1-C1' dihedral. Requires res0 and res1
        try:
            phi = self.computeDihedral(C0, N1, CA1, C1)
            if phi[1]==1:
                print('C0, N1 and CA1 are colinear in residue: ' + str(res1[0][5]))
        except:
            phi = [None, 0]

        #compute psi; the N1-Ca1-C1'-N2 dihedral. Requires res1 and res2
        try:      
            psi = self.computeDihedral(N1, CA1, C1, N2)
            if psi[1]==1:
                print('CA1, C1 and N2 are colinear in residue: ' + str(res1[0][5]))
        except:
            psi = [None, 0]
    
        #compute omega; the Ca1-C1-N2-CA2 dihedral. Requires res1 and res2
        try:
            omega = self.computeDihedral(CA1, C1, N2, CA2)
            if omega[1]==1:
                print('C1, N2 and CA2 are colinear in residue: ' + str(res1[0][5]))            
        except:
            omega = [None, 0]
        return (phi[0], psi[0], omega[0])
        
    # computes dihedral between four positions, returns angle, code.  If code ==1 then points are colinear. IF 0 calculation worked 
    def computeDihedral(self, P1, P2, P3, P4):
     
        # construct the in plane vectors V1 and U1 are both plane V and U respectively. UV is in both.
        U1 = P2-P1
        U1 = U1/np.linalg.norm(U1)
        UV = P3-P2
        UV = UV/np.linalg.norm(UV)
        V1 = P4-P3
        V1 = V1/np.linalg.norm(V1)
    
        # compute the dot product between the input vectors
        COSU = np.vdot(U1,UV)
        COSV = np.vdot(V1,UV)
    
        # set the default output value
        theta = 0.0
        colinear = 0
        # check for colinearity (COSU=1 or COSV=1)
        if (1.0-abs(COSU))>1E-6 and (1.0-abs(COSV))>1E-6:
            # Compute the normals to the planes
            N1 = np.cross(U1, UV)
            N1 = N1/np.linalg.norm(N1)
            N2 = np.cross(UV, V1)
            N2 = N2/np.linalg.norm(N2)
            
            # compute the binormal to N1 and UV.
            M1 = np.cross(UV,N1)
            
            # compute the components of N2 in the frame N1,UV,M1. The component of N2 along UV is always zero by definition.
            COSTHETA = np.vdot(N1,N2)
            SINTHETA = np.vdot(M1,N2)
    
            # shave off rounding and precision errors. COSTHETA should never be higher than 1 since N1 and N2 are normalised -doesn't matter too much in atan2 function.
            if COSTHETA>1.0:
                COSTHETA=1.0
    
            # Use the atan2 function which gives correct sign and a range between -180 and 180 in degrees,
            theta = np.degrees(np.arctan2(SINTHETA,COSTHETA))
        else:
            colinear = 1
    
        return theta, colinear
        
    