import numpy as np
import matplotlib.pyplot as plt
from signal_prep import generate_signal,spike_raster_plot_8

Fs = 22050.

# tone_1 = generate_signal(freq=1000,dBSPL=50.,duration=0.05,
#                        modulation_freq=0.,fs=Fs,ramp_duration=0.0025,plt=None,silence=True)
#
# for i,sample in enumerate(tone_1):
#     if abs(sample) > 1e-10:
#         tone_1 = tone_1[i:]
#         break

# tone_1 = tone_1[1104:]
#import SpiNNak-Ear output
input_directory  = '/home/rjames/Dropbox (The University of Manchester)/EarProject/Pattern_recognition/spike_trains/IC_spikes/'
output_directory = './'


n_fibres = '1000'
test_stimulus = '13.5_1_kHz'
duration = '75s'
intensity = '50dB'
test_file = 'spinnakear_'+test_stimulus+'_'+duration+'_'+intensity+'_'+n_fibres+'fibres'
an_data = np.load(input_directory+test_file+'.npz')
n_fibres_per_ihc = 10.

an_spikes = an_data['scaled_times']
n_ears= len(an_spikes.shape)

test_data_file = 'spinnakear_'+test_stimulus+'_38s_'+intensity+'_'+n_fibres+'fibres'
test_an_data = np.load(input_directory+test_data_file+'.npz')
test_an_spikes = test_an_data['scaled_times']
test_onset_times = test_an_data['onset_times']
#extract single test average spikes for the tone_1 stimulus




# spike_raster_plot_8(test_an_spikes, plt, 2., 1001, 0.001,
#                     title="test data activity")
# spike_raster_plot_8(an_spikes, plt, 2., 1001, 0.001,
#                     title="training data activity")
# plt.show()
#obtain stimulus onset times and durations
duration=150.
onset_times = an_data['onset_times']
n_stim = 0 #todo:sort for 2 ears
durations = [duration,duration]
# instantiate the auto correlation matrix for each stimulus based on the average of frequency bands across presentations
average_band_response=[]
for dur in durations:
    average_band_response.append(np.zeros((int((len(an_spikes)/n_fibres_per_ihc)),int(np.ceil(dur*0.001*Fs)))))

# for ear in range(n_ears):
#     test = [(1,2),(2,2)]
# n_reps = []
# for i,stimulus in enumerate(onset_times):
#     n_reps.append(0)
#     for time in stimulus:
#         n_reps[i]+=1
#         #find active
#         for j,neuron in enumerate(an_spikes):
#             neuron_response_times = np.asarray([int((t-time)*0.001*Fs) for t in neuron if t >= time and t < (time+duration)])
#             if len(neuron_response_times)>0:
#                 average_band_response[i][int(j/n_fibres_per_ihc),neuron_response_times]+=1
#         average_band_response[i] = average_band_response[i] / (n_reps[i] * np.ones(average_band_response[i].shape))
reverse_correlation_data = np.load(output_directory+test_file+'_reverse_correlation_data.npz')
# reverse_correlation_data = np.load(input_directory+test_file+'_reverse_correlation_data.npz')
average_band_response = reverse_correlation_data['average_band_response']
tone_1 = reverse_correlation_data['original_stimulus']

tone_1= np.append(tone_1,np.zeros(int(np.ceil(duration * 0.001 * Fs))-len(tone_1)))
#
tau_max = int(10 * 0.001 * Fs)
# for stim in average_band_response:
stim=average_band_response[1]
lag_matrix = np.zeros((len(stim)*tau_max, int(np.ceil(duration * 0.001 * Fs))))
for i,band in enumerate(stim):
    lag_matrix[(i*tau_max),:] = band[:]
    for j in range(1,tau_max):
        lag_matrix[(i*tau_max)+j,j:] = band[:-j]

auto_correlation = np.matmul(lag_matrix ,lag_matrix.transpose())
cross_correlation = np.matmul(lag_matrix,tone_1.transpose())
mapping_function = np.matmul(np.linalg.inv(auto_correlation),cross_correlation)

np.savez_compressed(input_directory+test_file+'_reverse_correlation_data.npz',original_stimulus=tone_1,average_band_response=average_band_response,
                   lag_matrix=lag_matrix,auto_correlation=auto_correlation,cross_correlation=cross_correlation,mapping_function=mapping_function)

print