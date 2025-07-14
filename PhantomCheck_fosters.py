# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:18:03 2025

@author: xanmc

General Analysis of Megin phantom
Checks HPI integrity, SNR, and localization for the standard phantom, either _raw or _ave FIF files.

For other visualization and processing options, see:
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import mne
import mne.utils
import matplotlib 
import mat73
import scipy.io

# Run to load Phantom Data to save for MATLAB processing #
raw_dir = '~/RESEARCH/Data/ILABS/Phantom/231207'
avg_file = ''
#raw_file = 'phantom_32_200nam_20240708_raw.fif'
raw_file = 'phantom_1000nam_default_IASoff_raw.fif'
raw = mne.io.read_raw_fif(os.path.join(raw_dir, raw_file),
                          preload=False, allow_maxshield='no')

## try saving as matlab matrix to fix time indexing issues
# data = raw.get_data()
# scipy.io.savemat("phantom_1000nam_default_IASoff_data.mat",{"raw_data":data,"raw_times":raw.times})

## uncomment to save covariance for Foster's in matlab
# cov_raw = mne.compute_raw_covariance(raw)
# scipy.io.savemat("phantom_32_200nam_20240708_raw_cov.mat",{"covar":cov_raw.data})

################
### uncomment to use matlab matrix to create new raw with our processed data
# datarec_fos = mat73.loadmat('1000nam_default_IASoff_check.mat')
# data = datarec_fos["phi_0_check"]
# ##sss
# datarec_vsh = mat73.loadmat('phantom_1000nam_default_IASoff_mat_ns_sss.mat')
# data = datarec_vsh["data_rec_vsh"]
# # # ##test one processed set at a time
# raw = mne.io.RawArray(data, correct_raw.info)
# del data #for memory'
##################

### try SSS with mne-python
#raw = mne.preprocessing.maxwell_filter(raw_correct, verbose=True)

# First, the frequencies: where they should be.
chpi_freqs, ch_idx, chpi_codes = mne.chpi.get_chpi_info(info=raw.info)
print(f'cHPI coil frequencies extracted from raw: {chpi_freqs} Hz')

# ##########################
##### Extract and examine phantom events
# Find events based on standard phantom trigger #
# Event IDs are embedded in 'SYS201' channel when using the standard Megin
# phantom stimulation protocol. Otherwise, need a custom algorithm.
events = mne.find_events(raw, stim_channel='SYS201', shortest_event=0)

# Create the epochs, assuming bursts of 2 20Hz cycles separated by .25 seconds #
# The baseline substraction should automatically be applied.
reject_dict = dict(grad=5000e-12, mag=15000e-12)
epochs = mne.Epochs(raw, events, tmin=-.05, tmax=0.15, baseline=(None, 0),
                    reject = reject_dict, preload=True, verbose=False)
print(f'{len(epochs)} epochs created from {len(events)} events.')
del raw #save space
## Display the waveforms as butterfly plots #
scalings = dict(mag=1e-11, grad=2e-10)
epochs.plot(butterfly=False, scalings=scalings);

##################
#### calculate and examine waveform averages
# Define averages #
# This creates a list of evoked objects, one for each phantom dipole location #
evokeds = epochs.average(picks='data', by_event_type=True);

## Plot the evoked waveforms (one dipole at a time)  butterfly style #
## Large SNR for all dipoles
idx = 20
fig, axset = plt.subplots(2, 1, figsize=[8, 6])  # for two axes side-by-side
evokeds[idx].plot(axes=axset, show=False);
fig.suptitle(f'Dipole {evokeds[idx].comment}', x=0.96, y=0.98,
              horizontalalignment='right', fontweight='bold');
fig.show();

# Distill each evoked object down to one value for each channel (per dipole) #
# No need to fit whole waveform. So could use single time point, max peak-to-peak
# amplitude, RMS energy, bandpass power, etc. 
# Method 1 - manually choose a high SNR time point; retains polarity
tpeak = .0375   # based on above average plots
ipeak = evokeds[0].time_as_index(.0375)[0]   # a useful, but don't need it here

evokeds_amp = []
for evk in evokeds:
    evk_new = evk.copy().crop(tpeak, tpeak)
    evokeds_amp.append(evk_new)

# Some setup, specific to the Megin/Elekta phantom  #
pos_actual, ori_actual = mne.dipole.get_phantom_dipoles()
sphere = mne.make_sphere_model(r0=(0., 0., 0.), head_radius=0.08)
cov = mne.make_ad_hoc_cov(evokeds[0].info)


# Make the dipole fits - just for the amplitude version here #
# Assuming that "true" phantom locations match order of dipoles activated
a = len(evokeds_amp)
a1 = pos_actual.shape[0]
assert len(evokeds_amp) == pos_actual.shape[0]   # expecting all to have been run

pos = np.nan * np.ones(pos_actual.shape)     # #dipoles x 3
ori = np.nan * np.ones(pos_actual.shape)
gof = np.nan * np.ones(pos_actual.shape[0])
dipole_results = []

print('Processing dipole... ', end='')
for ii, evk in enumerate(evokeds_amp):   # here, just for amplitude version
    assert evk.comment == f'{ii + 1}', evk.comment
    print(f'{evk.comment} ', end='')
    
    dip = mne.fit_dipole(evk, cov, sphere, verbose=False)[0]
    pos[ii] = dip.pos[0]
    ori[ii] = dip.ori[0]
    gof[ii] = dip.gof[0]
    dipole_results.append(dip)
print('\n')

pos_err = 1000 * np.linalg.norm(pos - pos_actual, axis=1)
ori_err = np.rad2deg(
    np.arccos(np.abs(np.sum(ori * ori_actual, axis=1))))

print(f'  Position error: mean={np.mean(pos_err):.1f} mm, '
      f'max={np.max(pos_err):.1f} mm')
print(f'  Orientation error: mean={np.mean(ori_err):.1f}°, '
      f'max={np.max(ori_err):.1f}°')
print(f'  Goodness of fit: mean={np.mean(gof):.1f}%, max={np.max(gof):.1f}%')











