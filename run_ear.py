import numpy
import model_launch_framework
#from mcmc_examples.lighthouse.lighthouse_model import LightHouseModel
from scipy.io import wavfile
import pylab as plt
import math
from vision.spike_tools.vis.vis_tools import plot_output_spikes
from signal_prep import *

def test_filter(audio_data,b0,b1,b2,a0,a1,a2):
    past_input=numpy.zeros(2)
    past_concha=numpy.zeros(2)
    concha=numpy.zeros(len(audio_data))
    for i in range(441,len(audio_data)):
        if i>=1202:
            print ''
        concha[i]=(b0 * audio_data[i]
                  + b1 * audio_data[i-1]#past_input[0]
                  + b2 * audio_data[i-2]#past_input[1]
                  - a1 * concha[i-1]#past_concha[0]
                  - a2 * concha[i-2]#past_concha[1]
                     ) * a0

        #past_input[1] = past_input[0]
        #past_input[0] = audio_data[i]

        #past_concha[1] = past_concha[0]
        #past_concha[0] = concha[i]
    return concha


#audio_data=numpy.fromfile("./c_models/load_files/load1_1_4k",dtype='float32')
#audio_data=numpy.fromfile("./c_models/load_files/load1_1_6k_22k",dtype='float32')
#audio_data=numpy.fromfile("./c_models/load_files/load1_1_6k_44k",dtype='float32')
#audio_data=numpy.fromfile("./c_models/load_files/load1_1_chirp",dtype='float32')
#audio_data=numpy.fromfile("./c_models/load_files/load1_1kate_22k",dtype='float32')
#audio_data=numpy.fromfile("./c_models/load_files/load1_1vowels_22k",dtype='float32')

audio_data = generate_signal(freq=4000,dBSPL=80.,duration=0.5,modulation_freq=5.)
numpy.save('../Brainstem/audio_data.npy',audio_data)

#concha = test_filter(audio_data,0.1783,0,-0.1783,1,-1.3477,0.6433)

#numpy.savetxt("./concha.csv",concha, fmt="%e", delimiter=",")

#audio_data=numpy.fromfile("./c_models/load_files/load1_1",dtype='float32')
#pole_freqs=numpy.fromfile("./c_models/load_files/pole_freqs_125",dtype='float32')
pole_freqs = numpy.logspace(3.3,3.95,50)#64)
#pole_freqs=[457, 6900]
#pole_freqs=numpy.fromfile("./c_models/load_files/pole_freqs",dtype='float32')
#pole_freqs=numpy.empty(25)#TODO:discover why this fails at 800+#
#pole_freqs.fill(8000.0)
#pole_freqs=numpy.empty(352)
#pole_freqs.fill(6900.0)
#pole_freqs=[30, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000]
#pole_freqs=[4000,4000]
#plt.figure()
#plt.plot(audio_data[12000:15000])#[12000:20500])
#
# plt.figure()
#plt.plot(audio_data)
#plt.ylim(-0.005,0.005)
#plt.show()

#check audio data can be divided evenly into 100 sample segements
audio_data = audio_data[0:int(numpy.floor(len(audio_data)/100)*100)]
#plt.figure()
#plt.plot(audio_data)
#plt.show

#create framework of connected model vertices and run
samples = model_launch_framework.run_model(
    audio_data, n_chips=numpy.ceil(len(pole_freqs)/2),n_drnl=2,
    pole_freqs=pole_freqs,n_ihcan=5,fs=22050,resample_factor=1,num_macks=4)
numpy.save("./samples.npy",samples)

#samples=numpy.load("./samples.npy")

#convert to spike train
spike_trains=list()
spike_index=0
ihc_index=0


#obtain list of IHCAN outputs
ihc_output = [samples[x:x+2*int(numpy.floor(len(audio_data)/100))*100] for x in xrange(0, len(samples), 2*int(numpy.floor(len(audio_data)/100))*100)]
drnl=numpy.zeros((len(pole_freqs)*10,int(numpy.floor(len(audio_data)/100))*100))
for ihc in ihc_output:
    #obtain fibre response
    hsr = [ihc[x:x+100] for x in xrange(0,len(ihc),200)]
    #HSR = [item for sublist in hsr for item in sublist]
    HSR=numpy.concatenate(hsr,axis=0)
    drnl[-spike_index][:]=HSR
    #find non zero indicies for fibre and record in spike trains list
    idxs = numpy.nonzero(HSR)
    if len(idxs[0]) == 0:
        print "spike_index {} from ihc{} did not record in full, re-simulate this particular channel".format(spike_index,ihc_index)
    for i in idxs[0]:
        spike_trains.append((spike_index, i))

    #increment spike_index
    spike_index=spike_index+1

    #obtain fibre response
    lsr=[ihc[x:x+100] for x in xrange(100,len(ihc),200)]
    LSR = numpy.concatenate(lsr,axis=0)
    drnl[-spike_index][:]=LSR
    #LSR = [item for sublist in lsr for item in sublist]
    # find non zero indicies for fibre and record in spike trains list
    idxs = numpy.nonzero(LSR)
    if len(idxs[0])==0:
        print "spike_index {} from ihc{} did not record in full, re-simulate this particular channel".format(spike_index,ihc_index)
    for i in idxs[0]:
        spike_trains.append((spike_index, i))

    # increment spike_index
    spike_index = spike_index + 1

    # increment ihc_index
    ihc_index=ihc_index + 1

#spike_trains=numpy.load("/home/rjames/Dropbox (The University of Manchester)/EarProject/spike_trains_kate_a_10kfib.npy")
duration = len(audio_data)/22050.0
spike_times = [spike_time for (neuron_id, spike_time) in spike_trains]
scale_factor = duration/numpy.max(spike_times)
scaled_times = [spike_time * scale_factor for spike_time in spike_times]
spike_ids = [neuron_id for (neuron_id, spike_time) in spike_trains]
spike_ids[:] = [neuron_id + 1 for neuron_id in spike_ids]

##plot results
plt.figure()
plt.plot(scaled_times, spike_ids, '.', markersize=1,
                 markerfacecolor='black', markeredgecolor='none',
                 markeredgewidth=0)

max_id = numpy.max(spike_ids)
size = len(spike_ids)
num_ticks=10.0
test=spike_ids[::int(numpy.floor((size)/num_ticks)-1)]
pole_freqs_size = len(pole_freqs)

ticks=[]
for i in range(int(num_ticks)+1):
    idx = int(numpy.floor((i / (num_ticks)) * (pole_freqs_size-1)))
    ticks.append(str(int(pole_freqs[idx])))

plt.yticks(test, ticks)

plt.ylim(numpy.min(spike_ids),numpy.max(spike_ids))
plt.xlim(0,numpy.max(scaled_times))
plt.xlabel('time (s)')
plt.ylabel('AN fibre best frequency (Hz)')
#plot_output_spikes(spike_trains,plotter=plt,markersize=1,color='black')'''
plt.show()

# Save the results
numpy.save("/home/rjames/Dropbox (The University of Manchester)/EarProject/spike_trains_6k_640fib_50dB.npy", spike_trains)
#numpy.savetxt("./results.csv",drnl, fmt="%e", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/results.csv", samples, fmt="%f", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/complete.txt",[1],fmt="%f")