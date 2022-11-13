#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 22:02:01 2022

@author: ms
"""

#%% Import modules
from PyQt5.QtWidgets import QApplication, QMainWindow, QTextEdit, QAction, QFileDialog
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import spikeinterface as si
import spikeextractors as se
import spiketoolkit as st
import spikesorters as ss
import spikecomparison as sc
import spikewidgets as sw
import neo
import scipy
#%% PyQt5
app = QApplication([])
#%% Variables
start_timestamp = 0
#%% Load data - Neo.io
# path = "/media/ms/2022_DATADRIVE_MS_2/ShuttleDrive/CheetahDataRegroup/221019_pcb3_2_1"
# print("\n\n"+path+"\n\n")
# # reader=neo.rawio.NeuralynxRawIO(path)
# recording=se.NeuralynxRecordingExtractor(dirname=path+'/')
# # reader.parse_header()
# recording.set_channel_locations([[0,602],[1,601],[1,602],[301,0],
#                                  [301,1],[302,0],[302,1],[301,301],
#                                  [301,302],[302,301],[0,0],[302,302],
#                                  [301,601],[301,602],[302,601],[302,602],
#                                  [601,0],[601,1],[602,0],[602,1],
#                                  [601,301],[0,1],[601,302],[602,301],
#                                  [602,302],[1,0],[1,1],[0,301],
#                                  [0,302],[1,301],[1,302],[0,601]])
# recording.set_channel_groups([2,2,2,3,
#                               3,3,3,4,
#                               4,4,0,4,
#                               5,5,5,5,
#                               6,6,6,6,
#                               7,0,7,7,
#                               7,0,0,1,
#                               1,1,1,2])
# # recording_prb = recording.load_probe_file(path+'/pcb5/1/probe1.prb')
# fs = recording.get_sampling_frequency()
# num_chan = recording.get_num_channels()
# channel_ids = recording.get_channel_ids()
# print('Channel ids:', channel_ids)
# print('Sampling frequency:', fs)
# print('Number of channels:', num_chan)
# # sos=scipy.signal.butter(20, [1,2000], btype='bandpass', fs=32000, output='sos')
# # recording_butter=scipy.signal.sosfilt(sos, recording)
# recording_ss=st.preprocessing.bandpass_filter(recording,freq_min=300,freq_max=8000,filter_type='butter')
# recording_cs=st.preprocessing.bandpass_filter(recording,freq_min=10,freq_max=8000,filter_type='butter')
#%% Load data
# path = QFileDialog.getExistingDirectory(directory='/media/ms/DATADRIVE_MS_1/[Exp]_ShuttleDrive/[Data]_ShuttleDrive_CheetahData',
#                                         caption='Select data directory')
path = "/media/ms/2022_DATADRIVE_MS_2/ShuttleDrive/CheetahData/2022-10-19_14-13-43/pcb5/1"
print("\n\n"+path+"\n\n")
recording = se.MdaRecordingExtractor(path)
channel_ids = recording.get_channel_ids()
fs = recording.get_sampling_frequency()
num_chan = recording.get_num_channels()
print('Channel ids:', channel_ids)
print('Sampling frequency:', fs)
print('Number of channels:', num_chan)
recording_prb = recording.load_probe_file(path+'/probe.prb')
recording_ss = st.preprocessing.bandpass_filter(recording,
                                                freq_min=300,
                                                freq_max=2000,
                                                filter_type='butter'
                                                )
recording_cs = st.preprocessing.bandpass_filter(recording,
                                                freq_min=1,
                                                freq_max=2000,
                                                filter_type='butter'
                                                )
# recording_ss = st.preprocessing.common_reference(recording_f_ss,
#                                                  reference='median'
#                                                  )
# recording_cs = st.preprocessing.common_reference(recording_f_cs,
#                                                   reference='median'
#                                                   )
#%% Mountainsort4_SS
sorting_SS = ss.run_mountainsort4(recording_ss,
                                  detect_sign=-1,
                                  # detect_threshold=3,
                                  filter=False,
                                  verbose=True,
                                  # num_workers=5,
                                  grouping_property='group'
                                  )
#%% Mountainsort4_CS
sorting_CS = ss.run_mountainsort4(recording_cs,
                                  detect_sign=1,
                                  # detect_threshold=1,
                                  filter=False,
                                  verbose=True,
                                  # num_workers=5,
                                  grouping_property='group'
                                  )
#%% re_units.xlsx
appended_data = []
appended_data = pd.DataFrame(appended_data)
for x in sorting_SS.get_unit_ids():
    print(x)
    a = sorting_SS.get_unit_spike_train(unit_id=x)
    a = a.tolist()
    a = pd.DataFrame(a,
                     columns=['SS_'+'{0:02d}'.format(x)]
                     )
    appended_data = pd.concat([appended_data, a],
                              ignore_index=False,
                              axis=1
                              )
for x in sorting_CS.get_unit_ids():
    print(x)
    a = sorting_CS.get_unit_spike_train(unit_id=x)
    a = a.tolist()
    a = pd.DataFrame(a,
                      columns=['CS_'+'{0:02d}'.format(x)]
                      )
    appended_data = pd.concat([appended_data, a],
                              ignore_index=False,
                              axis=1
                              )
appended_data /= fs
appended_data += start_timestamp
# appended_data.to_csv(path+'/units.csv')
appended_data.to_excel(path+'/units.xlsx')
#%% Post-processing
x = [recording_ss, recording_cs]
y = [sorting_SS, sorting_CS]
z = ['SS', 'CS']
for i, j, k in zip(x, y, z):
    templates = st.postprocessing.get_unit_templates(i,
                                                     j,
                                                     max_spikes_per_unit=200,
                                                     save_as_property=True,
                                                     verbose=True
                                                     )               
    max_chan = st.postprocessing.get_unit_max_channels(i,
                                                       j,
                                                       save_as_property=True,
                                                       verbose = True
                                                       )
    st.postprocessing.compute_unit_template_features(i,
                                                     j,
                                                     max_channels_per_waveforms=1
                                                     )
    se.MdaSortingExtractor.write_sorting(sorting=j,
                                         save_path=path+'/firings_'+k+'.mda'
                                         )
    st.postprocessing.export_to_phy(i,
                                    j,
                                    output_folder=path+'/phy_MS4_'+k
                                    )