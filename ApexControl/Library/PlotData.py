#!/usr/bin/env python3
'''
Created on 10 Mar 2020

@author: Chris
'''
import json
import csv
import numpy as np
import matplotlib.pyplot as plt

class plotData():
    
    def __init__(self, filename, dataType='EC', configFile=None):
        if configFile:
            self.settings = self.loadJson(configFile)['settings']
        else:
            self.settings = "UseDataSettings"
        
        if dataType=='EC':
            self.plotAll_IOT(self.loadECsv(filename))
        
        elif dataType=='IOT':
            self.plotAll_IOT(self.loadJson(filename))

        elif dataType=='PW':
            self.plotAll_PW(self.loadJson(filename))
        
        elif dataType=='Spect':
            self.plotSpectrum(self.loadJson(filename))
        
        elif dataType=='SLED':
            self.plotAll_SLED(self.loadJson(filename))

        elif dataType=='MultiSpec':
            self.plot_MultiSpectrum(self.loadJson(filename))

        elif dataType=='Formanator':
            self.plotAll_Formanator(self.loadJson(filename))

        elif dataType=='Bzometer':
            self.plotAll_Bzometer(self.loadJson(filename))
        
        elif dataType=='BzomEChem':
            self.plotAll_Bzometer(self.loadJson(filename))

        elif dataType=='Ramachandran':
            self.plotRamachandran(self.loadJson(filename))
        else:
            print("Unknown data type: ", dataType)

    def loadECsv(self, filename):
        with open( filename ,"r") as csvfile:
            data_csv = csv.reader(csvfile, delimiter=',')
            data_dict={}
            
            for row in data_csv:
                if len(row)==4:
                    try:
                        data_dict[ int( row[0] ) ][ 't' ].append( float( row[1] ) )
                        data_dict[ int( row[0] ) ][ 'v' ].append( float( row[2] ) )
                        data_dict[ int( row[0] ) ][ 'i' ].append( float( row[3] ) )
                    except KeyError:
                        data_dict[ int( row[0] ) ] = { 't': [ float( row[1] ) ], 'v': [ float( row[2] ) ], 'i': [ float( row[3] ) ] }

                if len(row)==3:
                    try:
                        data_dict[ 0 ][ 't' ].append( float( row[1] ) )
                        data_dict[ 0 ][ 'v' ].append( float( row[2] ) )
                        data_dict[ 0 ][ 'i' ].append( float( row[3] ) )
                    except KeyError:
                        data_dict[ 0 ] = {'t': [ float( row[1] ) ], 'v': [ float( row[2] ) ], 'i': [ float( row[3] ) ] }
                        
        return data_dict

    def loadJson(self, filename):
        """
        loads a test file from the test number 

        Args:
            testNum integer for the file
 
        """
        with open( filename ) as f:
            return json.load(f)

    def saveJson(self, filename, data):
        """
        Saves the data 

        Args:
            filename - the filename of the file
            data - the data to store in the file in json format

        """
        with open( filename, "w" ) as f:
            json.dump( data, f, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o.__dict__)


    # plots a ramachandran plot from the torsion
    def plotRamachandran(self, data):
        
        if self.settings=="UseDataSettings":
            settings = data['settings']
        else:
            settings = self.settings
        
        resInfo = data['resInfo']
        plotDict = settings['plotDict']
        plt.figure()
        
        for chain in resInfo:
            for res in resInfo[chain]['resList']:
                torsion = resInfo[chain][str(res)]['torsion']
                resName = resInfo[chain][str(res)]['atomList'][0][3]
                try:
                    plotString = plotDict[resName]
                except KeyError:
                    plotString = plotDict['DEF']
                if not torsion[0]==None and not torsion[1]==None:
                    plt.plot(torsion[0], torsion[1], plotString)
            
        plt.xlabel('Phi (degrees)')
        plt.ylabel('Psi (degrees)')
        if settings['Grid']==1:
            plt.grid()
        plt.xlim(settings['phiMin'], settings['phiMax'])
        plt.ylim(settings['psiMin'], settings['psiMax'])
        plt.title(settings['Title'])
        plt.show()

    def plotAll_ECsv(self, data_csv):
        for chan, data in data_csv.items():
            plt.figure(chan)
            plt.subplot(211)
            plt.plot(data['t'],data['v'])
            plt.ylabel('potential (V)')
            plt.grid('on')
            plt.title('volt,curr vs time, channel = {0}'.format(chan))
            
            plt.subplot(212)
            plt.plot(data['t'],data['i'])
            plt.ylabel('current (uA)')
            plt.xlabel('time (sec)')
            plt.grid('on')
            
            plt.figure(int(chan) + len(data_csv))
            plt.plot(data['v'],data['i'])
            plt.xlabel('potential (V)')
            plt.ylabel('current (uA)')
            plt.title('curr vs volt, channel = {0}'.format(chan))
            plt.grid('on')
        
            plt.show()

    def plotAll_IOT(self, data_dict):
        for chan, data in data_dict.items():
            plt.figure(chan)
            plt.subplot(211)
            plt.plot(data['t'],data['v'])
            plt.ylabel('potential (V)')
            plt.grid('on')
            plt.title('volt,curr vs time, channel = {0}'.format(chan))
            
            plt.subplot(212)
            plt.plot(data['t'],data['i'])
            plt.ylabel('current (uA)')
            plt.xlabel('time (sec)')
            plt.grid('on')
            
            plt.figure(int(chan) + len(data_dict))
            plt.plot(data['v'],data['i'])
            plt.xlabel('potential (V)')
            plt.ylabel('current (uA)')
            plt.title('curr vs volt, channel = {0}'.format(chan))
            plt.grid('on')
        
            plt.show()
            
            
    def plotSpectrum(self, testData):
        plt.figure()
        plt.plot(testData['wavelengths'], testData['spectrum'])      
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Intensity (A.U.)')
        plt.show()

    def plotAll_PW(self, testData):
        testNum = testData['testNum']

        for channel in range(1, len(testData['ECData']) + 1):
        
            ''' function takes the output of a set of a single test run and plots the information in an intelligible way.''' 
            t_EC = np.array(testData['ECData'][str(channel)]['t']) - testData['universalStartTime']
            i_EC = np.array(testData['ECData'][str(channel)]['i'])
            v_EC = np.array(testData['ECData'][str(channel)]['v'])
    
            wavelengths = np.array(testData['spectralData']['wavelengths'])
    
            t_spectra = np.array([ int(key) for key in testData['spectralData'] if key not in ['wavelengths'] ])
            
            spectra = np.array([ testData['spectralData'][str(t)] for t in t_spectra   ] )                 
            
        
            plt.figure(testNum)
            plt.subplot(311)
            plt.plot(t_EC, v_EC)
            plt.ylabel('potential (V)')
            plt.grid('on')
            plt.title('volt,curr vs time, test number = {0}, period = {1}, duty cycle ={2}, start time= {3}'.format(testNum, testData['pulsePeriod'], testData['dutyCycle'], testData['universalStartTime']))
        
            plt.subplot(312)
            plt.plot(t_EC, i_EC)
            plt.ylabel('current (uA)')
            plt.grid('on')
            plt.xlabel('time (sec)')
            
            plt.subplot(313)
            plt.plot(v_EC, i_EC)
            plt.xlabel('potential (V)')
            plt.ylabel('current (uA)')
            plt.grid('on')
    
            plt.figure(testNum + len(testData['ECData']))
            plt.imshow(spectra, aspect='auto', extent=[min(wavelengths), max(wavelengths), max(t_spectra) - testData['universalStartTime'], min(t_spectra) - testData['universalStartTime']])      
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('time (s)')
            
            plt.show()

    def plotAll_SLED(self, testData):
        testNum = testData['testNum']

        ''' function takes the output of a set of a single test run and plots the information in an intelligible way.''' 
        t = np.array(testData['MUSBData']['t']) - testData['universalStartTime']
        b = np.array(testData['MUSBData']['b'])
    
        wavelengths = np.array(testData['spectralData']['wavelengths'])
    
        t_spectra = np.array([ int(key) for key in testData['spectralData'] if key not in ['wavelengths'] ])
            
        spectra = np.array([ testData['spectralData'][str(t)] for t in t_spectra   ] )                 
            
        
        plt.figure(testNum)
        plt.plot(t, b)
        plt.ylabel('Brightness (0 to 100)')
        plt.xlabel('Time (s)')
        plt.grid('on')
        plt.title('brightness vs time, test number = {0}, period = {1}, amplitude ={2}, start time= {3}, phase offset= {4}'.format(testNum, testData['Period'], testData['Amplitude'], testData['universalStartTime'], testData['Phase']))
    
        plt.figure(testNum + len(testData['MUSBData']))
        plt.imshow(spectra, aspect='auto', extent=[min(wavelengths), max(wavelengths), max(t_spectra) - testData['universalStartTime'], min(t_spectra) - testData['universalStartTime']])      
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('time (s)')
        
        plt.show()


    def getTestString(self, dataDict):
        outStr = ""
        for key in dataDict:
            if key not in ['universalStopTime', 'universalStartTime','testNum']:
                outStr +=key
                outStr += "=" + str(dataDict[key]) + ", "
        
        return outStr + '\n'

    def plotEchem_BZometer(self, data):
        # plots out the data that is present in the given dictionary
        try:
            testString = self.getTestString(data['testDict'])
            testNum = data['testDict']['testNum']
        
            try:
                numChannels  = len(data['ECData'])

                for channel in range(1, numChannels + 1):
                    t_EC = np.array(data['ECData'][str(channel)]['t'])
                    i_EC = np.array(data['ECData'][str(channel)]['i'])
                    v_EC = np.array(data['ECData'][str(channel)]['v'])

                    plt.figure(channel)
                    plt.subplot(311)
                    plt.plot(t_EC, v_EC, '+')
                    plt.ylabel('potential (V)')
                    plt.xlabel('time (secs)')
                    plt.grid('on')
                    plt.title('volt, curr vs time, ' + testString + ' channel= ' + str(channel) + ', test= ' + str(testNum))

                    plt.subplot(312)
                    plt.plot(t_EC, i_EC, '+')
                    plt.ylabel('current (uA)')
                    plt.grid('on')
                    plt.xlabel('time (secs)')
                    
                    plt.subplot(313)
                    plt.plot(v_EC, i_EC, '+')
                    plt.xlabel('potential (V)')
                    plt.ylabel('current (uA)')
                    plt.grid('on')
                    plt.tight_layout()
                    plt.show()
            except KeyError as e:
                print("PlotData: testNum: ", data['testDict']['testNum'], " ECData Not Present: ", e)
            except TypeError as e:
                print("PlotData:  testNum: ", data['testDict']['testNum'], " ECData not populated: ", e)
        except KeyError:
            print("PlotData: Test Dictionary not present")
    
    def plotAll_Bzometer(self, data):
        # plots out the data that is present in the given dictionary
        try:
            testString = self.getTestString(data['testDict'])
            testNum = data['testDict']['testNum']
        
            try:
                numChannels  = len(data['ECData'])

                for channel in range(1, numChannels + 1):
                    t_EC = np.array(data['ECData'][str(channel)]['t'])
                    i_EC = np.array(data['ECData'][str(channel)]['i'])
                    v_EC = np.array(data['ECData'][str(channel)]['v'])

                    plt.figure(channel)
                    plt.subplot(311)
                    plt.plot(t_EC, v_EC, '+')
                    plt.ylabel('potential (V)')
                    plt.xlabel('time (secs)')
                    plt.grid('on')
                    plt.title('volt, curr vs time, ' + testString + ' channel= ' + str(channel) + ', test= ' + str(testNum))

                    plt.subplot(312)
                    plt.plot(t_EC, i_EC, '+')
                    plt.ylabel('current (uA)')
                    plt.grid('on')
                    plt.xlabel('time (secs)')
                    
                    plt.subplot(313)
                    plt.plot(v_EC, i_EC, '+')
                    plt.xlabel('potential (V)')
                    plt.ylabel('current (uA)')
                    plt.grid('on')
                    plt.tight_layout()
                    plt.show()
            except KeyError as e:
                print("PlotData: testNum: ", data['testDict']['testNum'], " ECData Not Present: ", e)
            except TypeError as e:
                print("PlotData:  testNum: ", data['testDict']['testNum'], " ECData not populated: ", e)
                
            try:
                wavelengths = np.array(data['spectralData']['wavelengths'])

                t_spectra = np.array(data['spectralData']['t'])
                t_spectra = np.floor(t_spectra).astype(int)
                spectra = np.array(data['spectralData']['spectra'])
                
                plt.figure()
                plt.imshow(spectra, aspect='auto', extent=[min(wavelengths), max(wavelengths), max(t_spectra) - data['testDict']['universalStartTime'], min(t_spectra) - data['testDict']['universalStartTime']])      
                plt.xlabel('Wavelength (nm)')
                plt.ylabel('time (s)')
                
                plt.show()
            except KeyError:
                print("PlotData:  testNum: ", data['testDict']['testNum'], " Spectral Data Not Present")
        except KeyError:
            print("PlotData: Test Dictionary not present")        


    def plot_MultiSpectrum(self, data):
        for spectrometerName in data:
            if spectrometerName not in ['universalStartTime']:
                specData = data[spectrometerName]
                wavelengths = np.array(specData['wavelengths'])
                t_spectra = np.array(specData['t'])
                t_spectra = np.floor(t_spectra).astype(int)
                spectra = np.array(specData['spectra'])
                
                plt.figure()
                plt.imshow(spectra, aspect='auto', extent=[min(wavelengths), max(wavelengths), max(t_spectra) - data['universalStartTime'], min(t_spectra) - data['universalStartTime']])      
                plt.xlabel('Wavelength (nm)')
                plt.ylabel('time (s)')
                plt.title('Spectrometer: ' + spectrometerName)
            
        plt.show()
    
    
    def plotAll_Formanator(self, tuningData):
        
        ''' function takes the output of a formanator PID and plots the information in an intelligible way.'''
        # First Generate Time Series
        t_keys = list(tuningData.keys())
        t_series = [ int(timeStamp) - int(t_keys[0]) for timeStamp in t_keys]
        spectra = [ np.array(tuningData[timeStamp]['spectrum']) for timeStamp in t_keys ]
        maxSpectra = [ max(spectrum) for spectrum in spectra]
        # pumpRate = [ float(tuningData[timeStamp]['pumpRate']) for timeStamp in t_keys ]
        # errVal = [ float(tuningData[timeStamp]['errVal']) for timeStamp in t_keys ]
        # errSum = [ float(tuningData[timeStamp]['errSum']) for timeStamp in t_keys ]
        wavelengthsDict = self.loadJson("CCSDefaultWavelengths.txt")
        wavelengths = wavelengthsDict['wavelengths']  
        configData = self.loadJson("config.json")
        setPoints = [float(configData['settings']['setPoint']) for _ in t_series]
        #waterRate = [np.abs(float(configData['settings']['pump1Rate'])) for _ in t_series]
   

        from matplotlib import rcParams
        rcParams['font.sans-serif'] = "Arial"
        # Then, "ALWAYS use sans-serif fonts"
        rcParams['font.family'] = "sans-serif"

        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12

        plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)

        fig = plt.figure(figsize=(3.6,3.6 * 3/4))
        plt.scatter(t_series, maxSpectra, marker='x', label='Peak Emission Intensity')
        plt.scatter(t_series, setPoints, marker='+', label='Set Point')
        plt.axvline(x=270, label="Small Manual Dilution", color='tomato')
        plt.axvline(x=301, label="Large Manual Dilution", color='red')
        plt.axvline(x=385, label="Recovery", color='green')
        plt.gca().grid(False)
        plt.ylabel('Peak Emission Intensity (arbitray unit)')
        plt.xlabel('Time (s)')
        plt.ylim((0, 0.5))
        #plt.xlim((0, 200))
        plt.title('Peak Emission Intensity vs Time')
        plt.legend()
        handles, labels = plt.gca().get_legend_handles_labels()
        order = [3, 4, 0, 1, 2]
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
        fig.savefig("PeakEmission.svg", format='svg', transparent=True)
        
        #fig = plt.figure()
        #plt.scatter(t_series, pumpRate, marker='x')
        #plt.title('Pump Rate vs Time')
        #plt.ylabel('pumpRate (steps per second)')
        #plt.xlabel('Time (s)')
        #plt.legend()
        #plt.grid('on')
        #fig.savefig("PumpRate.svg", format='svg')
        
        fig = plt.figure(figsize=(3.6,3.6 * 3/4))
        plt.imshow(spectra, aspect='auto', extent=[min(wavelengths), max(wavelengths), max(t_series), min(t_series)])      
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Time (s)')
        plt.title('Spectra vs Time')
        fig.savefig("Waterfall.svg", format='svg', transparent=True)
        
        plt.show()

if __name__ == '__main__':
    import sys
    filename = sys.argv[1]
    dataType = sys.argv[2]
    print(sys.argv)
    try:
        configFile = sys.argv[3]
    except IndexError:
        configFile = None
    
    print("Plotting filename: ", filename, " of data type: ", dataType, " with config file: ", configFile)
    plotObj = plotData(filename, dataType=dataType, configFile=configFile)
    print("Done")