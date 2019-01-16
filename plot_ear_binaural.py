import numpy as np
import model_launch_framework_binaural
from scipy.io import wavfile
import pylab as plt
import math
from signal_prep import *
from scipy.io import savemat, loadmat

# Save the results
results_directory = '/home/rjames/Dropbox (The University of Manchester)/EarProject/Pattern_recognition/spike_trains/IC_spikes'
results_file = "/spinnakear_timit_3s_30dB_30000fibres.npz"
# results_file = "/spinnakear_kate_a_31s_60dB_3000fibres.npz"
# results_file = "/spinnakear_a_21s_60dB_3000fibres.npz"
# results_file = "/spinnakear_yes_21s_60dB_3000fibres.npz"

results_data = np.load(results_directory+results_file)
binaural_scaled_times = results_data['scaled_times']
binaural_audio_data = results_data['audio_data']
Fs = results_data['Fs']

ear_list = ['left','right']
for ear,scaled_times in enumerate(binaural_scaled_times):
    duration = binaural_audio_data[ear].size/Fs
    spike_raster_plot_8(scaled_times, plt, duration, len(scaled_times) + 1, 0.001,
                        title="AN activity " + ear_list[ear] + " ear",markersize=1,subplots=(2,1,ear+1)
                        )#,filepath=results_directory)


plt.figure("input audio")
for channel in binaural_audio_data:
    t = np.linspace(0,duration,len(channel))
    plt.plot(t,channel,alpha=0.5)
    plt.xlabel("time (s)")
    plt.xlim((0,duration))
plt.legend(['left','right'])

plt.show()