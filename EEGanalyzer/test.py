from warnings import filters
import pyedflib, os
import pandas as pd
import numpy as np
from progress.bar import Bar
from scipy import signal
from setup import config

class EDF:

    filtersactive = 0

    def __init__(self, filepath):

        self.openFile(filepath)
        self.generate_signal_df()
        self.generate_annotation_df()
        self.find_trial_indices()
        self.initialize_filters()
        self.perform_operation('average')
        self.images = ProcessImages(self)

        EXCEL_BAR = Bar('Generating Spreadsheets', max=len(self.summaries))
        excel_filename = f'{config.new_csv_path}/{self.name} Operation Summary.xlsx'
        with pd.ExcelWriter(excel_filename) as writer:
            for key, value in self.summaries.items():
                value.to_excel(writer, sheet_name=f'Condition {key}')
                EXCEL_BAR.next()
            EXCEL_BAR.finish()
    
    def openFile(self, filepath):

        self.filedata = pyedflib.EdfReader(filepath)
        self.channels = self.filedata.getSignalLabels()
        self.num_channels = self.filedata.signals_in_file
        self.name = self.filedata.getPatientName()
        self.animal_number = self.name.split('-')[-1]

    def generate_annotation_df(self, annotation='TTL: Rise'):

        annotation_df = pd.DataFrame.from_records(self.filedata.readAnnotations()).transpose().reset_index(drop=True)
        annotation_df = annotation_df.drop(1, axis=1) #Remove duration
        annotation_df.columns = ['Time from Start', 'Annotation']
        annotation_df = annotation_df.loc[(annotation_df['Annotation'] == annotation)]
        annotation_df = annotation_df.drop_duplicates(subset='Time from Start', keep='first').reset_index(drop=True)
        annotation_df['Time from Start'] = annotation_df['Time from Start'].astype(float).round(4)

        self.label_pulses(annotation_df)
        self.annotations = annotation_df

    def generate_human_labels(self, num_conditions, num_trials):

        conditions = []
        trials = []
        for condition in np.arange(1,(num_conditions+1),1):
            for trial in np.arange(1,(num_trials+1),1):
                trials.append(trial)
                conditions.append(condition)
        
        return(conditions, trials)

    def find_trial_indices(self):

        trial_start_indices = []
        trial_end_indices = []

        TS_BAR = Bar('Labelling Indices', max=len(self.annotations))
        for trial_start in iter(self.annotations['Time from Start']):
            all_trials_after = self.waveforms.loc[(self.waveforms['Time from Start'] >= trial_start) & (self.waveforms['Time from Start'] < trial_start + 0.1)].index
            trial_start_index = int(all_trials_after[0])
            trial_start_indices.append(trial_start_index)
            trial_end_indices.append(int(trial_start_index + 667))
            TS_BAR.next()
        TS_BAR.finish()
        
        self.annotations['Start'] = trial_start_indices
        self.annotations['Stop'] = trial_end_indices
        
    def generate_signal_df(self):

        signal_dataframe = pd.DataFrame(columns=self.channels)
        for channel_id in range(0,self.num_channels):
            signal_dataframe[self.filedata.getLabel(channel_id)] = self.filedata.readSignal(chn=channel_id)
        signal_dataframe['Time from Start'] = self.label_waveform_time()
        
        self.waveforms = signal_dataframe

    def label_pulses(self, annotation_df):

        conditions, trials = self.generate_human_labels(15, 100)

        del conditions[annotation_df.shape[0]:]
        del trials[annotation_df.shape[0]:]

        annotation_df['Condition'] = conditions
        annotation_df['Trial'] = trials

        return(annotation_df)

    def label_waveform_time(self):

        end_time_seconds = self.filedata.getFileDuration()
        timestep = (self.filedata.getFileDuration() / self.filedata.getNSamples())
        if isinstance(timestep, np.ndarray):
            timestep = float(timestep[0])
        time_labels = np.arange(0,end_time_seconds,timestep)
        time_labels = np.around(time_labels,4)
        return(time_labels)

    def perform_operation(self, mode):
        
        condition_summaries = dict({})

        num_conditions = len(self.annotations['Condition'].unique())

        if(mode == 'average'):
            for condition in np.arange(1,(num_conditions+1),1):
                condition_summaries[condition] = (self.average_trials(condition))
        
        self.summaries = condition_summaries

    def initialize_filters(self):

        config.get_filters()
        self.bandpass = bool(config.bp_filter_init)
        self.notch = bool(config.notch_filter_init)

        if self.bandpass:
            self.sos = signal.butter(int(config.bp_quality_factor, [config.bp_low, config.bp_high], 'bandpass', fs=config.sample_freq, output='sos'))
            self.filtersactive += 1
        
        if self.notch:
            self.b, self.a = signal.iirnotch(config.notch_target, config.notch_quality_factor, config.sample_freq)
            self.filtersactive += 1

    def filter_data(self, raw_data=False, init=False):

        bandpass = self.bandpass
        notch = self.notch

        if notch and bandpass:
            notch_data = pd.Series(signal.filtfilt(self.b, self.a, raw_data))
            filtered_data = pd.Series(signal.sosfilt(self.sos, notch_data))
        
        if notch and not bandpass:
            filtered_data = pd.Series(signal.filtfilt(self.b, self.a, raw_data))
        
        if bandpass and not notch:
            filtered_data = pd.Series(signal.sosfilt(self.sos, raw_data))
        
        return(filtered_data)

    def average_trials(self, condition): #Called on each condition it is given
        
        #Get Annotation File, Waveforms, and Channel List
        annotations = self.annotations
        waveforms = self.waveforms
        channels = self.channels

        #Find all trials within the specified condition
        find_condition = (self.annotations['Condition'] == condition)
        trial_times = self.annotations.loc[find_condition].reset_index(drop=True)

        #Construct a Dataframe for Each Channel
        channel_dataframes = {}
        for channel in np.arange(0,len(channels), 1):
            channel_dataframes[channel] = pd.DataFrame()

        #Add/ Filter Raw Data for Each Trial
        for trial_id in np.arange(0, len(trial_times),1):

            row = trial_times.iloc[trial_id]
            start_index = int((row['Start']).astype(float))
            stop_index = int((row['Stop']).astype(float))

            #(by channel)
            for channel_id in np.arange(0, len(channels), 1):
                
                raw_data = self.waveforms.iloc[start_index:stop_index][channels[channel_id]].reset_index(drop=True)

                #Filter
                data = raw_data
                if self.filtersactive == 0:
                    data = raw_data
                else:
                    data = self.filter_data(raw_data)

                channel_data = channel_dataframes[channel_id]
                channel_data = channel_data.append(data, ignore_index=True)
                channel_data.index = np.arange(1, len(channel_data) + 1, 1)
                channel_dataframes[channel_id] = channel_data

        #Construct a summary of the data analysis
        summary_dict = {}
        
        for channel_id, dataframe in channel_dataframes.items():

            waveform = dataframe.mean()
            dataframe.loc['Average'] = waveform

            channel_dataframes[channel_id] = dataframe
            summary_dict[self.channels[channel_id]] = waveform
        
        summary_dataframe = pd.DataFrame.from_dict(summary_dict, orient='index')

        #Finally, return the condition summary
        return(summary_dataframe)

class ProcessImages():

    def __init__(self, analysis):
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
        except KeyError:
            raise KeyError('Not Configured for Image Output!')

        summary_dataframe = analysis.summaries
        outputer = []
        for key, value in summary_dataframe.items():
            outputer.append(key)
        print(outputer)
        pass

def main():
    filepath='C:/Users/Joseph/Documents/Data Analysis/EEG Data/Comprehensive/EEGanalyzer/EEGanalyzer/Data Files/2021-12-13/Whack.edf'
    test_edf = EDF(filepath)
    pass

if __name__ == '__main__':
    main()