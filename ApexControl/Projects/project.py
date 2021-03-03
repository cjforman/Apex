'''
Created on 1 Jul 2020

@author: Chris

Sets up the basic settings file handling of a project.

'''
import json
from FileIO import FileIO

class Project(FileIO):
    '''
    classdocs
    '''
    def __init__(self, settings):
        
        # copy any settings over
        self.settings = settings
        
        # TODO: create instance of all apex objects used in the project, 
        # which all take params dictionary as a constructor value
        #e.g. self.PS = PS.PotentiostatController(settings)

                
    @classmethod
    def from_config(cls, config):
        """
        Obtains the necessary information from a configuration setup.

        Args:
            cls (BZometer): The instantiating class.

            config (Dict): Dictionary containing the configuration data.

        """
        settings = config['settings']
        return cls( settings ) 

    @classmethod
    def from_configfile( cls, configfile ):
        """
        Obtains the configuration data from a configuration file.

        Args:
            cls (BZometer): The instantiating class.

            configfile (File): The configuration file.

        """
        with open( configfile ) as f:
            return cls.from_config( json.load( f ) )

