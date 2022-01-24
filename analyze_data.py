## VEP ANALYSIS TOOL V2.0 
# UNIVERSITY OF COLORADO ANSCHUTZ MEDICAL CAMPUS IN VIVO NEUROPHYSIOLOGY CORE
##
## PRIMARY DATA ANALYSIS MODULE (analyze_data.py)
##

# IMPORT REQUIRED FILES

from calendar import c
import wave
import matplotlib
from setup import fileHandler as fh, config, makeCSV, fileSystemInit, config, initFilter

from os import listdir, system
from os.path import isfile, join
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import signal
import re

from progress.bar import Bar, FillingCirclesBar

class Fileset:

    def __init__(self, experiment_files, run_files_dataframe):

        afn, afp, efn, efp = ['', '', '', '']

        friendly_name = re.split('\ |\#',experiment_files['Friendly Name'])
        self.experiment_date = friendly_name[0]
        self.animal_number = friendly_name[-1]

        for file in experiment_files['Files']:
            file = str(file).strip()

            # Find filetype by checking original dataframe
            file_info = run_files_dataframe[(run_files_dataframe['Filename'] == file)]
            file_type = str(file_info.iloc[0]['Type']).lower()
            file_path = str(file_info.iloc[0]['Filepath'])

            # Define annotation and export filenames
            if file_type == 'annotations':
                self.afn = file
                self.afp = file_path
            elif file_type == 'export':
                self.efn = file
                self.efp = file_path

    def loadfiles(self, afp, efp):

        FILELOAD = Bar('Generating Dataframes', max=2)

        #Read in Annotations File
        #Check filetype and select reader to populate dataframe
        at = pd.read_csv(afp, sep='\t', header=6)
        at.columns = ['Number', 'Start Time', 'End Time', 'Time from Start', 'Channel', 'Annotation']
        at = at[at.index < 4510]
        at['Number'] = at['Number'].astype(float)
        at['Number'] = at['Number'].astype(int)
        at['Start Time'] = pd.to_datetime(at['Start Time'], errors='coerce')
        at['End Time'] = pd.to_datetime(at['End Time'], errors='coerce')
        at['Time from Start'] = at['Time from Start'].astype(float).round(3)
        at = at.sort_values(by=['Time from Start'])

        FILELOAD.next()

        #Read in export file
        ts = pd.read_csv(efp, sep='\t', header=6)
        ts.columns = ['Date', 'Time', 'Time Stamp', 'Time from Start', 'EEG1', 'EEG2', 'EMG', 'TTL']
        ts['Date'] = pd.to_datetime(ts['Date'])
        ts['Time'] = pd.to_datetime(ts.Time, errors='coerce')
        ts['Time'] = ts.Time.dt.strftime('%H:%M:%S')
        ts['Time Stamp'] = ts['Time Stamp'].astype(float)
        ts['Time from Start'] = ts['Time from Start'].astype(float).round(3)
        ts['EEG1'] = ts['EEG1'].astype(float)
        ts['EEG2'] = ts['EEG2'].astype(float)
        ts['EMG'] = ts['EMG'].astype(float)

        self.at = at
        self.ts = ts

        FILELOAD.next()
        FILELOAD.finish()

    def genFiles(self):

        req_files = []

        try:
            req_folders = list(config.config['output-folders'])
        except:
            req_folders = ['Images', 'Spreadsheets']

        self.newpath = {}

        for folder in req_folders:
            filepath = join(config.output_path,self.experiment_date,self.animal_number,folder)
            req_files.append(filepath)
            self.newpath[folder] = filepath
        
        fileSystemInit(req_files)

    def getTTL(self, at, ts):
        
        trials = []
        conditions = []

        try:
            number_conditions = len(config.config['conditions'].keys())
            number_trials = int(config.config['trials-per-condition'])
        except:
            raise KeyError('Config Not Working')

        CONLAB_BAR = Bar('Labelling Trials', max=(number_conditions * number_trials))
        for x in range(1,(number_conditions+1)):
            for y in range(1,(number_trials + 1)):
                trials.append(y)
                conditions.append(x)
                CONLAB_BAR.next()
        CONLAB_BAR.finish()

        condition_iter = iter(conditions)
        trial_iter = iter(trials)

        TTL_times = at.loc[((at['Channel'] == "ALL") | (at['Channel'] == "EEG1")) & (at['Annotation'] == "TTL: Rise")].reset_index(drop=True)
        TTL_times['Condition'], TTL_times['Trial'], TTL_times['Start'], TTL_times['Stop'], TTL_times['Datapoints'] = [0, 0, 0, 0, 0]
        
        try:
            TTL_times['Condition'] = TTL_times['Condition'].apply(lambda x: next(condition_iter))
        except StopIteration:
            pass

        try:
            TTL_times['Trial'] = TTL_times['Trial'].apply(lambda x: next(trial_iter))
        except StopIteration:
            pass

        TTL_times = TTL_times.reset_index(drop=True)

        def gather_index_files(row, id, iterator):
        
            i = next(iterator)
            start_time = float(TTL_times['Time from Start'].iloc[i])

            try:
                end_time = float(TTL_times['Time from Start'].iloc[i+1])
            except:
                end_time = (start_time + 1)
            
            trial_index = ts.loc[(ts['Time from Start'] >= start_time) & (ts['Time from Start'] <= end_time)].index
            
            try:

                if (id == 'Start'):
                    return(trial_index[0])
                elif (id == 'Stop'):
                    return(trial_index[-1])
                elif (id == 'Datapoints'):
                    return((trial_index[-1] - trial_index[0]))

            except IndexError:
                raise IndexError('No Trials Were Identified')


        GI_BAR = Bar('Applying TTL Labels', max=3)
        iterator = iter(list(range(0,int(10e6))))
        TTL_times['Start'] = TTL_times['Start'].apply(gather_index_files, args=['Start', iterator])
        GI_BAR.next()
        iterator = iter(list(range(0,int(10e6))))
        TTL_times['Stop'] = TTL_times['Stop'].apply(gather_index_files, args=['Stop', iterator])
        GI_BAR.next()
        iterator = iter(list(range(0,int(10e6))))
        TTL_times['Datapoints'] = TTL_times['Datapoints'].apply(gather_index_files, args=['Datapoints', iterator])
        GI_BAR.next()
        GI_BAR.finish()

        self.TTL_times = TTL_times
        makeCSV(TTL_times, f'TTL_{self.efn}.csv', destination=self.newpath['Spreadsheets'])

    def Initialize(self):
    
        self.loadfiles(self.afp, self.efp)
        self.genFiles()
        self.getTTL(self.at, self.ts)
        run_analysis = Analysis(self)
        run_analysis.analyze_files()

class Analysis:

    def __init__(self, fileset):

        self.fileset = fileset

        try:
            self.sample_freq = config.config['sample-freq']
            self.numtrials = config.config['trials-per-condition']
            self.conditions = config.config['conditions']
            self.duration = config.config['trial-duration']
            self.channels = config.config['channel-names']
            self.fft_range = config.config['fft-range']
            self.generate_images = config.config['generate-images']

            self.frames_per_ms = float(self.sample_freq / 1000)
            self.numFrames = round(float(self.frames_per_ms * self.duration))
            self.numChannels = len(self.channels)
            self.numconditions = len(self.conditions)
            
        except KeyError:
            raise KeyError('Data Details Configured Improperly')
        
        self.filters_active = 0
        self.NotchFilterInit()
        self.BandpassFilterInit()
    
    def NotchFilterInit(self):
        active, target, quality_factor, sample_freq = initFilter('notch')

        if active:
            self.b, self.a = signal.iirnotch(target, quality_factor, sample_freq)
            self.filters_active += 1
            self.notch_filter_active = True

    def BandpassFilterInit(self):
        active, [bp_high, bp_low], quality_factor, sample_freq = initFilter('bandpass')

        if active:
            self.sos = signal.butter(int(quality_factor), [bp_low, bp_high], 'bandpass', fs=sample_freq, output='sos')
            self.filters_active += 1
            self.bandpass_filter_active = True

    def FilterData(self, raw_data):

        bandpass = bool()
        notch = bool()
        
        try:
            notch = self.notch_filter_active
        except AttributeError:
            notch = False
        
        try:
            bandpass = self.bandpass_filter_active
        except AttributeError:
            bandpass = False


        if (notch == True) and (bandpass == True):
            notch_data = pd.Series(signal.filtfilt(self.b, self.a, raw_data))
            filtered_data = pd.Series(signal.sosfilt(self.sos, notch_data))
        elif (notch == True) and (bandpass == False):
            filtered_data = pd.Series(signal.filtfilt(self.b, self.a, raw_data))
        elif (bandpass == True) and (notch == False):
            filtered_data = pd.Series(signal.sosfilt(self.sos, raw_data))
        else:
            raise KeyError('FilterData() thinks it should not see this')
        
        return(filtered_data)

    def spectrogram(self, waveform):

        sampling_rate = self.sample_freq
        T = (1 / sampling_rate)
        try:
            fft_range = [0, self.fft_range]
        except:
            raise ValueError('FFT Range Not Defined or not Integer')

        waveform_fft = np.fft.fft(waveform)
        N = waveform.size
        f = np.linspace(0, 1 / T, N)
        
        return(waveform_fft,N,f)

    def generateImages(self, summary, condition_number, condition_title):
        
        try:
            waveform_padding = config.config['waveform-padding']
            figure_width = config.config['figure-width']
            figure_height = config.config['figure-height']
            override_autoy = config.config['override-autoy']
            override_autoy_min = config.config['y-min']
            override_autoy_max = config.config['y-max']
            x_label = config.config['x-label']
            y_label = config.config['y-label']
            grid_pref = config.config['grid']
            fft_title = config.config['fft-title']
            fft_x_label = config.config['fft-x-label']
            fft_y_label = config.config['fft-y-label']
            extension = config.config['image-type']
            save_images = config.config['save-images']
            show_images = config.config['show-images']
            fft_override_autoy = config.config['fft-override-autoy']
            fft_magnitude_max = config.config['fft-magnitude-max']

        except:
            raise KeyError('Not Configured for Image Output')

        summary_dataframe = summary.transpose()

        figure_title = f'{self.fileset.experiment_date} Animal #{self.fileset.animal_number} Condition #{condition_number}: {condition_title}'

        fig, axs = plt.subplots(nrows=self.numChannels, ncols=2, sharex=False, sharey=False, figsize=(figure_width, figure_height))

        IMAGE_BAR = Bar('Generating Images', max=self.numChannels)
        for channel_id in np.arange(0,self.numChannels,1):

            channel_waveform = summary_dataframe.iloc[:,channel_id].values
            spectrum_waveform, N, f = self.spectrogram(channel_waveform)

            data_min = int(round(min(channel_waveform)))
            data_max = int(round(max(channel_waveform)))
            data_diff = (data_max - data_min)
            correction_factor = (data_diff / waveform_padding)
            channel_data_min = (data_min - correction_factor)
            channel_data_max = (data_max + correction_factor)

            if(override_autoy):
                channel_data_min = override_autoy_min
                channel_data_max = override_autoy_max

            index = np.arange(0,self.numFrames,1) #Set index to length specified in mS
            x_ticks = np.arange(0,self.duration,10)

            axs[channel_id, 0].set_title(f'{self.channels[(channel_id+1)]}')
            axs[channel_id, 0].plot(index, channel_waveform)
            axs[channel_id, 0].set_xlim([0,self.numFrames])
            axs[channel_id, 0].set_ylim([channel_data_min, channel_data_max])
            #axs[channel_id, 0].set_xticks(x_ticks)
            axs[channel_id, 0].set_xlabel(x_label)
            axs[channel_id, 0].set_ylabel(y_label)
            axs[channel_id, 0].grid(grid_pref)

            axs[channel_id, 1].set_title(fft_title)
            axs[channel_id, 1].plot(f[:N // 2], np.abs(spectrum_waveform)[:N // 2] * 1 / N)
            axs[channel_id, 1].set_xlim([0,self.fft_range])
            if(fft_override_autoy):
                axs[channel_id, 1].set_ylim([0,fft_magnitude_max])
            axs[channel_id, 1].set_xlabel(fft_x_label)
            axs[channel_id, 1].set_ylabel(fft_y_label)
            axs[channel_id, 1].grid(grid_pref)

            IMAGE_BAR.next()

        IMAGE_BAR.finish()
        fig.suptitle(figure_title)

        plt.tight_layout()

        if(show_images):
            plt.show()

        if(save_images):
            filename = f"{condition_number} {condition_title}"
            filepath = join(self.fileset.newpath['Images'], filename)
            plt.savefig(f'{filepath}.{extension}')
        
        return

    def analyze_files(self):
        
        condition = iter(self.conditions)
        for i in condition:
            
            TTL_times = self.fileset.TTL_times

            # Generate dataframes for raw data
            channel_dataframes = {}
            for channel in iter(self.channels):
                channel_dataframes[channel] = pd.DataFrame()

            # Find all TTL times matching current condition number
            findCondition = (TTL_times['Condition'] == i)
            
            #Locate within trial times the start and stop index of the condition
            trial_times = TTL_times.loc[findCondition].reset_index(drop=True)
            trial_iter = iter(trial_times.index)

            #For each trial found in trial times matching the condition
            CONDTRIAL_BAR = FillingCirclesBar(f'Condition {i}', max=self.numtrials)
            for trial_id in trial_iter:
                num_frames = self.numFrames #Framerate determined in __init__ and config
                row = trial_times.iloc[trial_id] #Find the row of the current trial to check TTL times for its indices
                start_ind = int((row['Start']).astype(float)) #Find the start index
                stop_ind = int(start_ind + num_frames) #Define stop index based on frequency and desired length

                #For all channels
                for channel_id in iter(self.channels):

                    channels = self.channels

                    ts = self.fileset.ts

                    raw_data = ts.iloc[start_ind:stop_ind][channels[channel_id]].reset_index(drop=True)

                    if (self.filters_active > 0):
                        data = self.FilterData(raw_data)
                    else:
                        data = raw_data
                    
                    channel_data = channel_dataframes[channel_id]
                    channel_data = channel_data.append(data, ignore_index=True)
                    channel_data.index = np.arange(1, len(channel_data) + 1)
                    channel_dataframes[channel_id] = channel_data

                CONDTRIAL_BAR.next()
            CONDTRIAL_BAR.finish()

            #Average Waveforms for Each Channel for the condition
            
            try:
                summary_dict = {}
                spect_dict = {}

                for channel_id, dataframe in channel_dataframes.items():

                    waveform = dataframe.mean()
                    dataframe.loc['Average'] = waveform

                    channel_dataframes[channel_id] = dataframe
                    #spect_dict[self.channels[channel_id]] = np.fft.fft(waveform)
                    summary_dict[self.channels[channel_id]] = waveform
                
            except:
                print(channel_dataframes)
                raise KeyError('No data coming out of the channels.')
            
            summary_page_dataframe = pd.DataFrame.from_dict(summary_dict, orient='index')
            #spectrum_page_dataframe = pd.DataFrame.from_dict(spect_dict, orient='index')

            EXCEL_BAR = Bar('Generating Spreadsheet', max=(len(channel_dataframes) + 1))
            excel_filepath = join(self.fileset.newpath['Spreadsheets'],f'{self.fileset.experiment_date}#{self.fileset.animal_number} Condition {i} {self.conditions[i]}.xlsx')
            with pd.ExcelWriter(excel_filepath) as writer:
                for key, value in channel_dataframes.items():
                    value.to_excel(writer, sheet_name=str(self.channels[key]))
                    EXCEL_BAR.next()
                summary_page_dataframe.to_excel(writer, sheet_name=f'Summary {self.fileset.animal_number} {self.fileset.experiment_date}')
                #spectrum_page_dataframe.to_excel(writer, sheet_name='Spectrogram (beta)')
                EXCEL_BAR.next()
            EXCEL_BAR.finish()
            
            if(self.generate_images):
                self.generateImages(summary_page_dataframe, i, self.conditions[i])
            else:
                pass

if __name__ == "__main__":

    raise KeyError('Run main.py!')