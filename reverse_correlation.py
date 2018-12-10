import numpy as np
import matplotlib.pyplot as plt
from signal_prep import generate_signal,spike_raster_plot_8
from scipy.io import loadmat, savemat , wavfile

def spikegram_generator(spikes,onset_times,duration,Fs,n_reps=None):
    if n_reps == None:
        n_reps = len(onset_times)
    spikegram = np.zeros((len(spikes)*n_reps,int(duration*0.001*Fs)))
    overlap_spikes = [[] for _ in range (len(spikes)*n_reps)]
    # for rep_index,time in enumerate(onset_times):
    for rep_index in range(n_reps):
        # offset_time = time - (2 * (tau_max/Fs)*1000.)
        offset_time = onset_times[rep_index] - (2 * (tau_max/Fs)*1000.)
        #find active neurons within time window
        for j,neuron in enumerate(spikes):
            neuron_response_times = np.asarray([int((t-offset_time)*0.001*Fs)
                                                for t in neuron if t >= offset_time and t < (offset_time+duration)])
            if neuron_response_times.size>0:
                #add firing times spikegram
                spikegram[j*n_reps+rep_index,neuron_response_times]+=1
                overlap_spikes[j*n_reps+rep_index] = neuron_response_times
    return spikegram,overlap_spikes

def average_band_response_generator(spikes,onset_times,tau_max,duration,Fs,n_fibres_per_ihc=10.):
    # average_band_response = np.zeros((int((len(spikes)/n_fibres_per_ihc)),int(np.ceil(duration*0.001*Fs))))
    # average_band_response = np.zeros((int((len(spikes)/n_fibres_per_ihc)),int(duration*0.001*Fs)))
    # n_reps=0
    # for time in onset_times:
    #     offset_time = time - (2 * (tau_max/Fs)*1000.)
    #     n_reps+=1
    #     #find active neurons within time window
    #     for j,neuron in enumerate(spikes):
    #         neuron_response_times = np.asarray([int((t-offset_time)*0.001*Fs)
    #                                             for t in neuron if t >= offset_time and t < (offset_time+duration)])
    #         if neuron_response_times.size>0:
    #             #add firing times to a count for all fibres in associated frequency band
    #             average_band_response[int(j/n_fibres_per_ihc),neuron_response_times]+=1
    # #average firing counts across number of presentations
    # # average_band_response /= (n_reps * n_fibres_per_ihc * np.ones(average_band_response.shape))
    # average_band_response /= (n_fibres_per_ihc * np.ones(average_band_response.shape))

    average_band_response = np.zeros((int((len(spikes)/n_fibres_per_ihc)),int(duration*0.001*Fs)))
    for j,neuron in enumerate(spikes):
        neuron_response_times = np.asarray([int(t*0.001*Fs) for t in neuron if t<duration])
        if neuron_response_times.size>0:
            #add firing times to a count for all fibres in associated frequency band
            average_band_response[int(j/n_fibres_per_ihc),neuron_response_times]+=1
    #average firing counts across number of presentations
    average_band_response /= (n_fibres_per_ihc * np.ones(average_band_response.shape))
    return  average_band_response

def lag_matrix_generator(stim,tau_max,duration,step=10):
    num_entry_rows = len(range(0,int(tau_max),step))
    lag_matrix = np.zeros((len(stim) * num_entry_rows, int(duration * 0.001 * Fs)))
    for i,band in enumerate(stim):
        lag_matrix[(i*num_entry_rows),:] = band[:]
        for k,j in enumerate(range(step,tau_max,step)):
            lag_matrix[(i*num_entry_rows)+k+1,j:] = band[:-j]
    return lag_matrix

def stimulus_reconstruction(stimulus_length,tau_max,test_stim,mapping_function,step=10):
    num_entry_rows = len(range(0,int(tau_max),step))
    S_t = np.zeros(int(stimulus_length))
    for t in range(int(tau_max), int(S_t.size)):
        for n in range(len(test_stim)):
            for g_index,tau in enumerate(range(0,int(tau_max),step)):
                S_t[t] += mapping_function[(n * num_entry_rows) + g_index] * test_stim[n][t - tau]
    return S_t


#import SpiNNak-Ear output
input_directory  = '/home/rjames/Dropbox (The University of Manchester)/EarProject/Pattern_recognition/spike_trains/IC_spikes/'
output_directory = './'

Fs = 22050.
n_fibres = '3000'
test_stimulus = 'yes'
duration = '21s'
intensity = '60dB'
test_file = 'spinnakear_'+test_stimulus+'_'+duration+'_'+intensity+'_'+n_fibres+'fibres'
an_data = np.load(input_directory+test_file+'.npz')
n_fibres_per_ihc = 10.
an_spikes = an_data['scaled_times']
onset_times = an_data['onset_times']
n_ears= len(an_spikes.shape)#todo:sort for 2 ears
test_data_file = 'spinnakear_'+test_stimulus+'_2s_'+intensity+'_'+n_fibres+'fibres'
test_an_data = np.load(input_directory+test_data_file+'.npz')
test_an_spikes = test_an_data['scaled_times']
test_onset_times = test_an_data['onset_times']
test_stimulus = test_an_data['audio_data'][0]

#extract single test average spikes for the tone_1 stimulus
tau_max_ms = 10
tau_max = int(tau_max_ms * 0.001 * Fs)
step = 10#int(1 * 0.001 * Fs)#int(tau_max/10.)
#todo: the full binaural sound file can be obtained from the spinnakear sim output

#trim stimulus so initial sample is onset time - 2*tau_max
stimulus = an_data['audio_data'][0]
onset_sample = int(((onset_times[0][0][0]*0.001) * Fs) - 2 * tau_max)
offset_sample = int(((onset_times[0][0][5]*0.001) * Fs) - 2 * tau_max)
stimulus = stimulus[onset_sample:offset_sample]
duration=(len(stimulus)/Fs)*1000. #ms

test_onset_sample = int(((test_onset_times[0][0][0]*0.001) * Fs) - 2 * tau_max)
test_offset_sample = int(((test_onset_times[0][0][1]*0.001) * Fs) - 2 * tau_max)
test_stimulus = stimulus[test_onset_sample:test_offset_sample]
test_duration = (len(test_stimulus)/Fs)*1000. #ms

# reverse_correlation_data = np.load(input_directory+test_file+'_reverse_correlation_data.npz')
# s_t = reverse_correlation_data['reconstructed_stimulus']
# wavfile.write(output_directory+test_file+'_reconstruction.wav',rate=Fs,data=s_t)

# plt.figure()
# plt.plot(s_t)
# plt.show()
# mapping_function=reverse_correlation_data['mapping_function']
# mapping_function = loadmat(input_directory+'g.mat')['g']
# average_band_response=reverse_correlation_data['average_band_response']
# average_band_response_test=reverse_correlation_data['average_band_response_test']
# mat_dict={}
# spikegram,overlap_spikes = spikegram_generator(an_spikes[0],onset_times[0][0],duration,Fs,n_reps=len(onset_times[0][0]))
# test_spikegram,test_overlap_spikes = spikegram_generator(test_an_spikes[0],test_onset_times[0][0],duration,Fs,n_reps=1)
# spike_raster_plot_8(overlap_spikes,plt,duration/1000.,len(overlap_spikes),scale_factor=1./Fs)
# spike_raster_plot_8(test_overlap_spikes,plt,duration/1000.,len(test_overlap_spikes),scale_factor=1./Fs)
# # spike_raster_plot_8(an_spikes[0],plt,duration/1000.,len(an_spikes[0]))
# plt.show()
# mat_dict['original_stimulus'] = stimulus
# mat_dict['spikegram'] = spikegram
# mat_dict['spikegram_test']= test_spikegram
# savemat('/home/rjames/av_response',mat_dict)
# print "n_reps={}".format(len(onset_times[0][0]))

# savemat('/home/rjames/av_response',reverse_correlation_data)
# plt.figure('s_t')
# plt.plot(s_t)
# plt.plot(tone_1)
# plt.show()

average_band_response = average_band_response_generator(an_spikes[0],onset_times[0][0],tau_max,duration,Fs,n_fibres_per_ihc=10)
average_band_response_test = average_band_response_generator(test_an_spikes[0],[test_onset_times[0][0][0]],tau_max,test_duration,Fs,n_fibres_per_ihc=10)
#
plt.figure()
plt.imshow(average_band_response,vmin=abs(average_band_response).min(), vmax=abs(average_band_response).max(), extent=[0, 1, 0, 1])
plt.colorbar()
plt.figure()
plt.imshow(average_band_response_test,vmin=abs(average_band_response_test).min(), vmax=abs(average_band_response_test).max(), extent=[0, 1, 0, 1])
plt.colorbar()
plt.show()
print "built average response matrices"
lag_matrix = lag_matrix_generator(average_band_response,tau_max,duration,step=step)
print "built lag matrix"
auto_correlation = np.matmul(lag_matrix ,lag_matrix.transpose())
print "completed auto correlation"
cross_correlation = np.matmul(lag_matrix,stimulus.transpose())
print "completed cross correlation"
mapping_function = np.linalg.solve(auto_correlation,cross_correlation)
print "generated mapping function"
# S_t = stimulus_reconstruction(duration*Fs*0.001,tau_max,average_band_response_test,mapping_function,step=step)
S_t = stimulus_reconstruction(test_duration*Fs*0.001,tau_max,average_band_response_test,mapping_function,step=step)
print "reconstructed stimulus"
plt.figure('s_t')
plt.plot(S_t)

plt.figure('stimulus')
plt.plot(test_stimulus)
# np.savez_compressed(input_directory+test_file+'_reverse_correlation_data.npz',original_stimulus=tone_1,average_band_response=average_band_response,
#                    lag_matrix=lag_matrix,auto_correlation=auto_correlation,cross_correlation=cross_correlation,mapping_function=mapping_function,
#                     reconstructed_stimulus = S_t)

# np.savez_compressed(input_directory+test_file+'_reverse_correlation_data.npz',original_stimulus=tone_1,
#                    mapping_function=mapping_function,
#                     reconstructed_stimulus = S_t)
wavfile.write(input_directory+test_file+'_reconstruction_lag{}ms_step{}.wav'.format(tau_max_ms,step),rate=Fs,data=S_t)
print "complete"
plt.show()

