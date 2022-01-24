## VEP ANALYSIS TOOL V2.0 
# UNIVERSITY OF COLORADO ANSCHUTZ MEDICAL CAMPUS IN VIVO NEUROPHYSIOLOGY CORE
##
## FRONT-END (main.py)
##

# IMPORT LIBRARIES FOR FRONT-END
from numpy import sort
import os
import re
import pandas as pd

import setup
from progress.bar import Bar
from progress.spinner import Spinner
from setup import fileHandler as fh, config, fileSystemInit, makeCSV
import analyze_data as analysis

class RunFolder:

    def __init__(self):

        ## Select folder to be analyzed [Run]
        if (config.config['constant-run-folder'] == False):
            run_folder_filepath = fh.folder_selection_popup()
        else:
            run_folder_filepath = config.config['run-folder-path']

        self.rff = run_folder_filepath
        self.subfolders = fh.findSubfolders(self.rff)
        self.runfiles = fh.FiltFiles(self.subfolders)
        self.populate()
    
    def populate(self):

        dataframe_elements = list(config.config['output-elements'])

        #Create dataframe from the dict
        run_files_dataframe = pd.DataFrame.from_dict(self.runfiles, orient='index').reset_index()
        run_files_dataframe.columns = ['Filename', 'Filepath'] #First 2 columns need to be set first as the dict keys

        #Add all user configured 'output_elements' to run_files dataframe
        for element in dataframe_elements:
            run_files_dataframe[str(element)] = str('')

        ## Parse the filename and populate the run files dataframe
        def sortRunfiles(row, operation, iterator):

            i = next(iterator)

            get_row = run_files_dataframe.iloc[i] #Gets the particular row from the dataframe

            filename = str(get_row['Filename']) #Pull these values from dataframe
            filepath = str(get_row['Filepath'])

            parsed_filename = re.split('-|_|,|\.|\:|\ ', filename) #Parses - _ , . : ' '

            filename_elements = list(config.config['filename-elements']) #Gets user filename elements

            fd = {} #Filename Dict

            new_iterator = iter(list(range(0,int(10e6))))
            #Parses for each element
            for element in filename_elements:
                j = next(new_iterator)
                fd[element] = parsed_filename[j]

            #Parses the time-data into strings - HARDCODED
            surgery_date = str(f"{fd['Surgery Year']}-{fd['Surgery Month']}-{fd['Surgery Day']}")
            experiment_date = str(f"{fd['Experiment Year']}-{fd['Experiment Month']}-{fd['Experiment Day']}")
            experiment_time = str(f"{fd['Experiment Hour']}:{fd['Experiment Minute']}:{fd['Experiment Second']}")

            #Parses Friendly Name
            friendly_name = f"{experiment_date} {experiment_time} Animal #{fd['Animal Number']}"

            #This stack handles which operation we are performing
            if(operation == 'Type'):
                return(str(fd['Type']).upper())
            elif(operation == 'Experiment Name'):
                return(fd['Experiment Name'])
            elif(operation == 'Experiment Date'):
                return(experiment_date)
            elif(operation == 'Experiment Time'):
                return(experiment_time)
            elif(operation == 'Animal Number'):
                return(fd['Animal Number'])
            elif(operation == 'Surgery Date'):
                return(surgery_date)
            elif(operation == 'Friendly Name'):
                return(friendly_name)
            else:
                return("Error")

        COLUMN_BAR = Bar('Sorting Run Files', max=len(dataframe_elements))
        for column in dataframe_elements:
            # Applies the function to populate the dataframe from user-configured dataframe_elements
            iterator = iter(list(range(0,int(10e6))))
            run_files_dataframe[column] = run_files_dataframe[column].apply(sortRunfiles,operation=column, iterator=iterator)
            COLUMN_BAR.next()
        COLUMN_BAR.finish()

        self.run_files_dataframe = run_files_dataframe
        ## Sort the dataframe to find pairs of files (export and annotations for a particular recording)
        unique_files = run_files_dataframe['Friendly Name'].unique() #Find unique friendly names
        run_sets = {}
        #For each friendly name
        UN_BAR = Bar('Matching Loose Socks', max=len(unique_files))
        for friendly_name in unique_files:
            #Pull all rows matching the friendly name
            file_set = run_files_dataframe.loc[(run_files_dataframe['Friendly Name'] == friendly_name)]
            #Only pull ones which have 2 matches
            if (len(file_set) == 2):
                file_names = list(file_set['Filename']) #List all filenames
                file_types = list(file_set['Type']) #List all filetypes

                #Filter out data where both files are the same type
                if(file_types[0] != file_types[1]):
                    #Add to run_sets
                    run_sets[friendly_name] = [file_names]
                else:
                    pass   
            else:
                pass
            UN_BAR.next()
        UN_BAR.finish()

        analysis_files_dataframe = pd.DataFrame.from_dict(run_sets, orient='index').reset_index()
        analysis_files_dataframe.columns = ['Friendly Name', 'Files']
        self.analysisfiles = analysis_files_dataframe
        makeCSV(self.analysisfiles, f'AnalysisFile.csv')

    def RunAnalysis(self):
        ANAL_SPINNER = Spinner('Analyzing Runfolder')
        for row, file_set in self.analysisfiles.iterrows():
            fileset = analysis.Fileset(file_set, self.run_files_dataframe)
            fileset.Initialize()
            ANAL_SPINNER.next()
        ANAL_SPINNER.finish()

if __name__ == "__main__":
    
    try:
        run_folder = RunFolder()
    except:
        raise KeyError('Something Went Wrong!')

    run_folder.RunAnalysis()