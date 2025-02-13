# Apex
Automated Parametric EXplorer—controlling hardware, data acquisition and analysis in closed loop systems

This is the PDB analysis code and the plot data routine. 

To install, clone the git repository somewhere on your system.  Make sure the Script directory is on your system path and the ApexControl, ApexControl/Library and ApexControl/Project directories are on your PythonPath for which ever python interpreter you are using. 

You can then navigate to a working directory containing your PDB file and create a "config.json" file which controls the parameters of the analysis.  There is a usage documentation included in the source code which gets printed when the code is missing a parameter in the config file. 

That config file can be called anything and here is an example, which is included in the directory for giggles.

{
    "settings": {
	"PDBFilename": "MASP2.pdb",
    "command": "resInfo",
	"atomTypes": ["ATOM"],
    "outputFileroot": "MASP2",
    "resInfoType": ["torsion"],
	"plotData": ["Ramachandran"],
    "verbose": 0,
	"heartBeat": 1000,
	"plotDict": {"PRO":"g+", "TYR": "bx", "DEF": "rx"},
	"Grid": 1,
	"phiMin": -180, 
	"phiMax":  180, 
	"psiMin": -180, 
	"psiMax":  180,
	"Title": "MASP2 Ramachandran Plot"
        }
}


To run the analysis script on windows powershell, when in the working directory containing the PDB file and the config.json file type in 

doPDBAnalyzer.ps1 .\config.json

and it should call your python interpreter passing the appropriate functions and data to the script for analysis. 

The analysis routine that is called is defined by "command" parameter. At present only the "resInfo" command is implemented wich loads and parses the PDB and creates a data structure (python dictionary) broken down by chains and residues.  It is designed to have various switches in the resInfoType, such as 'torsion', which are different kinds of analyses. The output from these analyses are added to the data structure against the residue or chains of interest. Right now all it does is add torsion information to the python dictionary at the residue level.  BUt we can expand that if we want to and switch different things on or off, including them easily in the config file construct.

When the analysis is finished the routine will save a json of the data structure. The config.json information with all the parameters is included in that saved data structure. So if you change stuff and give it new file names then the parameters used in the data are saved along side it. 

The PDB Proc will then call the plot data command which loads the combined saved json and config data and calls the list of plotting functions given in the plotData list. At the moment only Ramachandran is implemented. But we can add more kinds of plots. 

You have to save the output image manually from those windows but its easy to add commands to save the figures. 

If you like once you have analysed the data once and created a json of the analysis, you can call the plot data command directly from the command line and re plot the data without running the analysis again. At that point you can control the plotting parameters in a separate plotting config json, or by omitting the plotting_config.json filename you can force the use of the plotting parameters defined in the saved data, created at analysis time from the original config.json.  The plot command does not over write the saved data ever. 

PlotData.ps1 dataFile.json Ramachandran plotting_config.json

That about does it. 

We can talk about what new functions and analysis we can add using this framework.

You should try a short PDB and make sure that the torsion output is the same as the way it is when you measure it in pymol. It should be but there are definitions that might need adjusting to guarantee that. 

Thanks



