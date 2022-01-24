## VEP ANALYSIS TOOL V2.0 
# UNIVERSITY OF COLORADO ANSCHUTZ MEDICAL CAMPUS IN VIVO NEUROPHYSIOLOGY CORE
##
## SETUP MODULE (setup.py)
##

import pathlib, yaml, os
from tkinter.tix import MAIN
from random import randint #For loading and parsing configuration
from tkinter import Tk #For generating GUI elements
from tkinter.filedialog import askdirectory #For file pop-up

config = ''

MAIN_FILEPATH = pathlib.WindowsPath(__file__).parent.parent.resolve()
curr_dir_filepath = pathlib.WindowsPath(__file__).parent.resolve() #Get Current Directory (WINDOWS)
config_filepath = os.path.join(MAIN_FILEPATH,'docs', 'config.yaml') #Configure Filepath to Config

def fileSystemInit(folders):
    try:
        dirs_made = 0

        FOLDERS = list(folders)
        
        for filepath in FOLDERS:
            if not os.path.exists(filepath):
                try:
                    oldmask = os.umask(0)
                    os.makedirs(filepath, 0o755)
                    dirs_made += 1
                finally:
                    os.umask(oldmask)
            else:
                pass
    except:
        raise PermissionError('FOLDER CREATION FAILED MAY NOT BE PERMISSIONS THO SORRY')

class Config:

    global MAIN_FILEPATH
    global curr_dir_filepath
    global config_filepath

    def __init__(self):

        self.MAIN_FILEPATH = pathlib.WindowsPath(__file__).parent.parent.resolve()
        self.curr_dir_filepath = pathlib.WindowsPath(__file__).parent.resolve() #Get Current Directory (WINDOWS)
        self.config_filepath = os.path.join(MAIN_FILEPATH,'docs', 'config.yaml') #Configure Filepath to Config
        self.config = self.read_yaml()
        self.default_extension = '.tsv'

        self.define_filepaths()
        self.generateFiles()
    
    def read_yaml(self):
        try:
            with open(self.config_filepath, "r") as f:
                return yaml.safe_load(f)
        except:
            raise FileExistsError('Configuration File Unable to be Loaded')
    def define_filepaths(self):
        ## DEFINE OUTPUT FILEPATHS FROM CONFIG FILE OR USE DEFAULT
        if(self.config['custom-output-folder']) != False:
            self.output_path = os.path.join(self.config['output-folder-path'])
        else:
            self.output_path = os.path.join(curr_dir_filepath, 'Output')

        if(self.config['custom-csv-folder']) != False:
            self.new_csv_path = os.path.join(self.config['new-csv-path'], '/')
        else:
            self.new_csv_path = os.path.join(self.output_path,'CSV')

        if self.config['use-default-file-extension'] != False:
            self.default_extension = '.tsv'
        elif self.config['analysis-file-extension']:
            self.default_extension = self.config['analysis-file-extension']
        else:
            self.default_extension = '.tsv'

        self.default_run_folder = os.path.join(curr_dir_filepath, 'Run')
    def generateFiles(self):
        fileSystemInit([self.default_run_folder, self.output_path, self.new_csv_path])
    def get_filters(self):
        try:
            self.sample_freq = self.config['sample-freq']
            self.notch_filter_init = self.config['notch-filter']
            self.notch_target = self.config['notch-target']
            self.notch_quality_factor = self.config['notch-quality-factor']

            self.bp_filter_init = self.config['bandpass-filter']
            self.bp_high = self.config['bandpass-high']
            self.bp_low = self.config['bandpass-low']
            self.bp_quality_factor = self.config['bandpass-quality-factor']
        except:
            raise AttributeError('Filters Not Setup in Configuration!')

if(isinstance(config, Config)):
    pass
else:
    config = Config()

class fileHandler:

    def __init__(self):
        print("File Handler Initialized")

    def folder_selection_popup(self):
        Tk().withdraw()
        try:
            folder_name = askdirectory()
        except:
            raise KeyError('AHHHHHHHH')
        if not folder_name:
            raise NameError('File Selection Error!')
        else:
            return(folder_name)

    def findSubfolders(self, filepath):

        subfolders = {}

        for subfolder in os.scandir(filepath):

            if subfolder.is_dir():
                subfolders[subfolder.name] = subfolder.path

        return(subfolders)

    def findSubfiles(self, filepath, extension=''):

        subfiles = {}

        for subfile in os.scandir(filepath):
            if subfile.is_file():
                filename = str(subfile.name)
                if extension == '':
                    subfiles[subfile.name] = subfile.path
                elif filename.endswith(extension):
                    subfiles[subfile.name] = subfile.path
        return(subfiles)

    def FiltFiles(self, subfolders, extension=config.default_extension):

        filt_files = {}
        for key, filepath in subfolders.items():
            add_files = self.findSubfiles(filepath, extension)
            filt_files.update(add_files)
        
        return(filt_files)

def makeCSV(dataframe, filename='', extension='.csv', destination=config.new_csv_path):

    try:
        if not filename:
            num = randint([1,100])
            filename = 'Untitled_{num}{extension}'
        
        filepath = os.path.join(destination, filename)
        dataframe.to_csv(filepath)
    except FileNotFoundError:
        raise FileNotFoundError(f'Output Folder Not Found! {destination}')
    except PermissionError:
        raise KeyError('All Fi')
    except:
        raise KeyError('CSV Failed to Come Into Existence')

def initFilter(filter):
    
    if(filter == 'notch'):
        try:
            sample_freq = config.sample_freq
            notch_filter_init = config.notch_filter_init
            notch_target = config.notch_target
            notch_quality_factor = config.notch_quality_factor

            return(notch_filter_init, notch_target, notch_quality_factor, sample_freq)
        except:
            raise KeyError('Notch Filter Not Complete in Config.yaml')
    if(filter == 'bandpass'):
        try:
            sample_freq = config.sample_freq
            bp_filter_init = config.bp_filter_init
            bp_high = config.bp_high
            bp_low = config.bp_low
            bp_quality_factor = config.bp_quality_factor

            return(bp_filter_init, [bp_high, bp_low], bp_quality_factor, sample_freq)

        except:
            raise KeyError('Bandpass Filter Not Complete in Config.yaml')
    else:
        raise KeyError(f'Filter {filter} not instantiated')


if __name__ == '__main__':
    'Run Main.py!'
else:
    fileHandler = fileHandler()
    config = Config()