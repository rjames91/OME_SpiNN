import numpy as np

target_duration_ms = 10.*60.*1000.
results_directory = "/home/rjames/Dropbox (The University of Manchester)/EarProject/Pattern_recognition/spike_trains/IC_spikes"

an_spikes=np.load(results_directory+'/ic_spikes_asc_train_60s.npy')

#find duration
max_time = 0

for neuron in an_spikes:
    if neuron.size>0 and neuron.max() > max_time:
        max_time = neuron.max()

#calculate how many repetitions we need
n_repeats = np.ceil(target_duration_ms / max_time)

ascending_audio_train_spikes = []#np.asarray([[]for _ in range(an_spikes.size)])

for id,neuron in enumerate(an_spikes):
    # if id==0:
    #     ascending_audio_train_spikes=np.asarray([j*max_time+time for j in xrange(int(n_repeats)) for time in neuron])
    # else:
    #     b=np.asarray([j*max_time+time for j in xrange(int(n_repeats)) for time in neuron])
    #     ascending_audio_train_spikes=np.vstack((ascending_audio_train_spikes,b))
    ascending_audio_train_spikes.append([j*max_time+time for j in xrange(int(n_repeats)) for time in neuron])

np.savez(results_directory + '/long_spike_trains', ascending_audio_train_spikes=ascending_audio_train_spikes)