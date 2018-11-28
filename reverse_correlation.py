import numpy as np
import matplotlib.pyplot as plt
from signal_prep import generate_signal,spike_raster_plot_8
from scipy.io import loadmat

def average_band_response_generator(spikes,onset_times,tau_max,duration,Fs,n_fibres_per_ihc=10.):
    # average_band_response = np.zeros((int((len(spikes)/n_fibres_per_ihc)),int(np.ceil(duration*0.001*Fs))))
    average_band_response = np.zeros((int((len(spikes)/n_fibres_per_ihc)),int(duration*0.001*Fs)))
    n_reps=0
    for time in onset_times:
        offset_time = time - (2 * (tau_max/Fs)*1000.)
        n_reps+=1
        #find active neurons within time window
        for j,neuron in enumerate(spikes):
            neuron_response_times = np.asarray([int((t-offset_time)*0.001*Fs)
                                                for t in neuron if t >= offset_time and t < (offset_time+duration)])
            if neuron_response_times.size>0:
                #add firing times to a count for all fibres in associated frequency band
                average_band_response[int(j/n_fibres_per_ihc),neuron_response_times]+=1
    #average firing counts across number of presentations
    average_band_response /= (n_reps * n_fibres_per_ihc * np.ones(average_band_response.shape))
    return  average_band_response

def lag_matrix_generator(stim,tau_max,duration):
    # lag_matrix = np.zeros((len(stim) * tau_max, int(np.ceil(duration * 0.001 * Fs))))
    lag_matrix = np.zeros((len(stim) * tau_max, int(duration * 0.001 * Fs)))
    for i,band in enumerate(stim):
        lag_matrix[(i*tau_max),:] = band[:]
        for j in range(1,tau_max):
            lag_matrix[(i*tau_max)+j,j:] = band[:-j]
    return lag_matrix

def stimulus_reconstruction(stimulus_length,tau_max,test_stim,mapping_function):
    S_t = np.zeros(int(stimulus_length))
    for t in range(int(tau_max), int(S_t.size)):
        for n in range(len(test_stim)):
            for tau in range(int(tau_max)):
                S_t[t] += test_stim[n][t - tau] * mapping_function[(n * tau_max) + tau]
    return S_t


#import SpiNNak-Ear output
input_directory  = '/home/rjames/Dropbox (The University of Manchester)/EarProject/Pattern_recognition/spike_trains/IC_spikes/'
output_directory = './'

Fs = 22050.
n_fibres = '1000'
test_stimulus = '13.5_1_kHz'
duration = '75s'
intensity = '50dB'
test_file = 'spinnakear_'+test_stimulus+'_'+duration+'_'+intensity+'_'+n_fibres+'fibres'
an_data = np.load(input_directory+test_file+'.npz')
n_fibres_per_ihc = 10.
an_spikes = an_data['scaled_times']
onset_times = an_data['onset_times']
n_ears= len(an_spikes.shape)#todo:sort for 2 ears
test_data_file = 'spinnakear_'+test_stimulus+'_38s_'+intensity+'_'+n_fibres+'fibres'
test_an_data = np.load(input_directory+test_data_file+'.npz')
test_an_spikes = test_an_data['scaled_times']
test_onset_times = test_an_data['onset_times']
#extract single test average spikes for the tone_1 stimulus
tau_max = int(10 * 0.001 * Fs)
#todo: the full binaural sound file can be obtained from the spinnakear sim output
tone_1 = generate_signal(freq=1000,dBSPL=50.,duration=0.05,
                       modulation_freq=0.,fs=Fs,ramp_duration=0.0025,plt=None,silence=True)
#trim stimulus so initial sample is onset time - 2*tau_max
# tone_1 = tone_1[int((onset_times[1][0]*0.001*Fs)-2*tau_max):]
for i,sample in enumerate(tone_1):
    if abs(sample) > 1e-10:
        tone_1 = tone_1[int(i - 2 * tau_max):]
        break
duration=(len(tone_1)/Fs)*1000. #ms

reverse_correlation_data = np.load(input_directory+test_file+'_reverse_correlation_data.npz')
s_t = reverse_correlation_data['reconstructed_stimulus']
# mapping_function=reverse_correlation_data['mapping_function']
mapping_function = loadmat(input_directory+'g.mat')['g']
average_band_response=reverse_correlation_data['average_band_response']
average_band_response_test=reverse_correlation_data['average_band_response_test']

plt.figure('s_t')
plt.plot(s_t)
# plt.plot(tone_1)
# plt.show()

# average_band_response = average_band_response_generator(an_spikes,onset_times[1],tau_max,duration,Fs)
# average_band_response_test = average_band_response_generator(an_spikes,[onset_times[1][0]],tau_max,duration,Fs)
# lag_matrix = lag_matrix_generator(average_band_response,tau_max,duration)
# auto_correlation = np.matmul(lag_matrix ,lag_matrix.transpose())
# cross_correlation = np.matmul(lag_matrix,tone_1.transpose())
# # mapping_function = np.matmul(np.linalg.inv(auto_correlation),cross_correlation)
# mapping_function = np.linalg.lstsq(auto_correlation,cross_correlation)
# reverse_correlation_data = np.load(output_directory+test_file+'_reverse_correlation_data.npz')
# average_band_response = reverse_correlation_data['average_band_response']
# mapping_function = reverse_correlation_data['mapping_function']
# auto_correlation =reverse_correlation_data['auto_correlation']
# cross_correlation =reverse_correlation_data['cross_correlation']
# lag_matrix = reverse_correlation_data['lag_matrix']

S_t = stimulus_reconstruction(duration*Fs*0.001,tau_max,average_band_response_test,mapping_function)
plt.figure('s_t matlab')
plt.plot(S_t)
plt.show()
# np.savez_compressed(input_directory+test_file+'_reverse_correlation_data.npz',original_stimulus=tone_1,average_band_response=average_band_response,
#                    lag_matrix=lag_matrix,auto_correlation=auto_correlation,cross_correlation=cross_correlation,mapping_function=mapping_function,
#                     reconstructed_stimulus = S_t)

# np.savez_compressed(input_directory+test_file+'_reverse_correlation_data.npz',original_stimulus=tone_1,
#                    mapping_function=mapping_function,
#                     reconstructed_stimulus = S_t)

print
