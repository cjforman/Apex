'''
Created on 14 Oct 2020

@author: Chris

Handles all file IO functions for APEX.
The projects module extends this class, so any class which is a project automatically has fileIO capability. 
Some of the individual libraries that aren't projects also extend this class, such as plotData, so that all the fileIO functionality is contained here. 

'''
import json
import csv
import cv2 as cv
import numpy as np
import sys
from collections.abc import Mapping

class FileIO(object):
    '''
    classdocs
    '''


    def __init__(self, params):
        '''
        Constructor
        '''
        self.params = params
    
    
    def saveObject(self, objectId, configfile):
        """
        Saves the configuration data for the given object

        Args:
            configfile (File): The configuration filename

        """
        with open( configfile, "w" ) as f:
            json.dump( objectId, f, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o.__dict__)

    def saveData(self, filename, data, delimiter=' '):
        """
        Saves data. 
        If the data is a mapping object, saves the data as a json.
        otherwise assumes data is a list of parameter arrays and attempts to save as a csv   

        Args:
            filename - the filename of the file
            data - the data to store in the file
            delimiter defaults to a space for a csv type file
        """
        if isinstance(data, Mapping):
            with open( filename, "w" ) as f:
                json.dump( data, f, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o.__dict__)
        else:
            with open( filename, mode='w') as csvfile:
                filewriter = csv.writer(csvfile, delimiter=delimiter)
                for row in data:
                    filewriter.writerow(row)

        return

    def loadData(self, filename, fileType='csv', imMode=cv.IMREAD_GRAYSCALE):
        # loads data either as csv or json type as specified. Defaults to csv
        if fileType =='csv':
            data = []
            with open( filename, newline='') as csvfile:
                filereader = csv.reader(csvfile, delimiter=' ')
                for row in filereader:
                    data.append(row)

        # implmented        
        elif fileType=='json':
            with open( filename) as f:
                data = json.load(f)
                
        elif fileType=='tif':
            data = cv.imread(filename, imMode)
        
        elif fileType=='txt':
            data = self.readTextFile(filename)
            
        else:
            data = {}

        return data

    def readTextFile(self, filename):
        #read line data in from file
        try:
            vst = open(filename, 'r')
        except IOError as e:
            print("I/O error", e.errno, e.strerror)
            print("Unable to open input file: " + filename)
            sys.exit()
        lines = vst.readlines()
        vst.close()
        return lines

    def writeTextFile(self, lines, filename):
        #write line data to file
        try:
            vst = open(filename, 'w')
        except:
            raise Exception("Unable to open output file: " + filename)
        for line in lines:
            a=line
            vst.write(a)
        vst.close()
        return

