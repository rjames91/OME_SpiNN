import numpy
import model_launch_framework
from scipy.io import wavfile
import pylab as plt
import math
from signal_prep import *
from scipy.io import savemat, loadmat


Fs = 22050.#44100.#24000.#34000.#10000-0.#40000.#
seg_size = 96
bitfield = True#False#
profile = False#True
resample_factor = 1.
#audio_data=numpy.fromfile("./c_models/load_files/load1_1",dtype='float32')
#audio_data=numpy.fromfile("./c_models/load_files/load1_1_6k_22k",dtype='float32')
#audio_data=numpy.fromfile("./c_models/load_files/load1_1_6k_44k",dtype='float32')
#audio_data=numpy.fromfile("./c_models/load_files/load1_1_chirp",dtype='float32')
#audio_data=numpy.fromfile("./c_models/load_files/load1_1kate_22k",dtype='float32')
#audio_data=numpy.fromfile("./c_models/load_files/load1_1vowels_22k",dtype='float32')
#audio_data = wavfile.read('./vowels_44_1kHz.wav')
#audio_data = numpy.float32(audio_data[1])
#audio_data=numpy.fromfile("/home/rjames/Dropbox (The University of Manchester)/EarProject/map_sig",
#                          dtype='float64')

# audio_data = generate_signal(freq=1000,dBSPL=20.,duration=1.,
#                              modulation_freq=0.,fs=Fs,ramp_duration=0.0025,plt=None,silence=True,silence_duration=0.005)
#audio_data = numpy.concatenate((audio_data,audio_data,audio_data))

# matlab = loadmat("/home/rjames/Dropbox (The University of Manchester)/EarProject/MAP_BS/map_64.mat")
#
# audio_data = matlab["map_output"][0]
#
# mat  =loadmat("/home/rjames/Dropbox (The University of Manchester)/EarProject/MAP_BS/map_cdispl.mat")
# audio_data = mat["map_output"][0]

"""[fs, raw_audio]= wavfile.read('./no_samples/0bde966a_nohash_1.wav')
plt.figure()
plt.plot(raw_audio)
plt.show()
cut_audio= raw_audio[0:8000]
wavfile.write('./no_samples/no_edit10.wav',fs,cut_audio)"""

# audio_data = generate_signal(signal_type="sweep_tone", freq=[30, 8000], dBSPL=50., duration=0.5,
#                             modulation_freq=0., fs=Fs, ramp_duration=0.0025, plt=None, silence=True)

dBSPL = 20.
sweep_duration = 0.1

asc = generate_signal(signal_type="sweep_tone", freq=[30, 8000], dBSPL=dBSPL, duration=sweep_duration,
                            modulation_freq=0., fs=Fs, ramp_duration=0.0025, plt=None, silence=True,ascending=True)

des = generate_signal(signal_type="sweep_tone", freq=[30, 8000], dBSPL=dBSPL, duration=sweep_duration,
                            modulation_freq=0., fs=Fs, ramp_duration=0.0025, plt=None, silence=True,ascending=False)

a = generate_signal(signal_type='file',dBSPL=dBSPL,fs=Fs,ramp_duration=0.0025,silence=True,
                            file_name='./a.wav',plt=None)
i = generate_signal(signal_type='file',dBSPL=dBSPL,fs=Fs,ramp_duration=0.0025,silence=True,
                            file_name='./i.wav',plt=None)
u = generate_signal(signal_type='file',dBSPL=dBSPL,fs=Fs,ramp_duration=0.0025,silence=True,
                            file_name='./u.wav',plt=None)

sounds = [asc,des]#[asc,des,a,i,u]#

audio_data = []
required_total_time = 60.
onset_times = [[]for _ in range(len(sounds))]
#for _ in range(n_repeats):
while int(len(audio_data)/Fs) < required_total_time:
    rand_choice = numpy.random.randint(len(sounds))
    chosen_samples = sounds[rand_choice]
    onset_found = False
    for sample in chosen_samples:
        if not onset_found and abs(sample)>1e-10:
            onset_times[rand_choice].append(numpy.round(1000.*(len(audio_data)/Fs)))
            onset_found = True
        audio_data.append(sample)


#plt.show()
#audio_data = numpy.concatenate((audio_data,audio_data,audio_data))

#audio_data = generate_signal(signal_type='file',dBSPL=40.,fs=Fs,ramp_duration=0.01,silence=True,silence_duration=0.1,
#                             file_name='./10speakers_2numbers_5repeats.wav',plt=None)#10speakers_2numbers_5repeats.wav

# audio_data = generate_signal(signal_type='file',dBSPL=40.,fs=Fs,ramp_duration=0.01,silence=False,silence_duration=0.1,
#                             file_name='./yes_samples/yes_edit1.wav',plt=None)#10speakers_2numbers_5repeats.wav

#plt.show()
#numpy.save('../Brainstem/audio_data.npy',audio_data)
#concha = test_filter(audio_data,0.1783,0,-0.1783,1,-1.3477,0.6433)
#numpy.savetxt("./concha.csv",concha, fmt="%e", delimiter=",")
#audio_data=numpy.fromfile("./c_models/load_files/load1_1",dtype='float32')
#pole_freqs=numpy.fromfile("./c_models/load_files/pole_freqs_125",dtype='float32')

if Fs > 34000.:
    pole_freqs = numpy.logspace(1.477,4.25,100)#full hearing spectrum
    rt = False
else:
    pole_freqs = numpy.logspace(1.477,3.95,100)#up to 9k range
    rt= True
   # rt = False

#pole_freqs=numpy.empty(8)
#pole_freqs.fill(9000)
#pole_freqs[0:4]=1000
#pole_freqs[4:10]=9000
#pole_freqs=numpy.empty(16)
#pole_freqs.fill(4000.0)
#pole_freqs=[30, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000]
#pole_freqs=[1000,1000]
#plt.figure()
#plt.plot(audio_data[12000:15000])#[12000:20500])
#
# plt.figure()
#plt.plot(audio_data)
#plt.ylim(-0.005,0.005)
#plt.show()

#check audio data can be divided evenly into 100 sample segements
audio_data = audio_data[0:int(numpy.floor(len(audio_data)/seg_size)*seg_size)]
#plt.figure()
#plt.plot(audio_data)
#plt.show
last_non_zero = numpy.nonzero(audio_data)[0].max()
#print "audio sample {}:{:.100e}".format(last_non_zero,audio_data[last_non_zero].astype(numpy.float32))


# plt.figure()
# plt.plot(audio_data)
#plt.show()

#create framework of connected model vertices and run
samples = model_launch_framework.run_model(
    audio_data, n_chips=numpy.ceil(len(pole_freqs)/2),n_drnl=2,
    pole_freqs=pole_freqs,n_ihcan=5,fs=Fs,resample_factor=resample_factor,num_macks=4,
    bitfield=bitfield,rt=rt,profile=profile)
#numpy.save("./samples.npy",samples)

#samples=numpy.load("./samples.npy")

#convert to spike train
spike_trains=[]
spike_times = []
for _ in range(len(pole_freqs)*10):
    spike_times.append([])
spike_index=0
ihc_index=0


#obtain list of IHCAN outputs
ihc_output = [samples[x:x+2*int(numpy.floor(len(audio_data)/seg_size))*seg_size]
              for x in xrange(0, len(samples), 2*int(numpy.floor(len(audio_data)/seg_size))*seg_size)]
drnl=numpy.zeros((len(pole_freqs)*10,int(numpy.floor(len(audio_data)/seg_size))*seg_size),dtype=numpy.float32)
for ihc in ihc_output:
    #obtain fibre response
    hsr = [ihc[x:x+seg_size] for x in xrange(0,len(ihc),seg_size*2)]
    HSR=numpy.concatenate(hsr,axis=0)
    drnl[-spike_index][:]=HSR
    if bitfield:
        #find non zero indicies for fibre and record in spike trains list
        idxs = numpy.nonzero(HSR)
        #if len(idxs[0]) == 0:
            #print "spike_index {} from ihc{} did not record in full, re-simulate this particular channel".format(spike_index,ihc_index)
        # for i in idxs[0]:
        #     spike_trains.append((spike_index, i))
        spike_times[spike_index]=idxs[0]

    #increment spike_index
    spike_index=spike_index+1

    #obtain fibre response
    lsr=[ihc[x:x+seg_size] for x in xrange(seg_size,len(ihc),seg_size*2)]
    LSR = numpy.concatenate(lsr,axis=0)
    drnl[-spike_index][:]=LSR
    if bitfield:
        # find non zero indicies for fibre and record in spike trains list
        idxs = numpy.nonzero(LSR)
        #if len(idxs[0])==0:
            #print "spike_index {} from ihc{} did not record in full, re-simulate this particular channel".format(spike_index,ihc_index)
        # for i in idxs[0]:
        #     spike_trains.append((spike_index, i))
        spike_times[spike_index]=idxs[0]

    # increment spike_index
    spike_index = spike_index + 1
    # increment ihc_index
    ihc_index=ihc_index + 1

if bitfield:
    #spike_trains=numpy.load("/home/rjames/Dropbox (The University of Manchester)/EarProject/spike_trains_kate_a_10kfib.npy")
    #[spike_trains,scale_factor]=numpy.load("./spike_trains.npy")
    duration = len(audio_data)/Fs
    # spike_times = [spike_time for (neuron_id, spike_time) in spike_trains]
    max_time = 0
    for neuron in spike_times:
        if neuron.size > 0:
            if neuron.max()>max_time:
                max_time = neuron.max()
    scale_factor = duration/max_time
    scaled_times = [1000*spike_time * scale_factor for spike_time in spike_times]
    spike_raster_plot_8(scaled_times, plt, duration, len(pole_freqs)*10 + 1, 0.001, title="pre pop activity")
    plt.show()

    # spike_ids = [neuron_id for (neuron_id, spike_time) in spike_trains]
    # spike_ids[:] = [neuron_id + 1 for neuron_id in spike_ids]
    #
    # ##plot results
    # plt.figure()
    # plt.plot(scaled_times, spike_ids, '.', markersize=1,
    #                  markerfacecolor='black', markeredgecolor='none',
    #                  markeredgewidth=0)
    #
    # max_id = numpy.max(spike_ids)
    # size = len(spike_ids)
    # num_ticks=10.0
    # test=spike_ids[::int(numpy.floor((size)/num_ticks)-1)]
    # pole_freqs_size = len(pole_freqs)
    #
    # ticks=[]
    # for i in range(int(num_ticks)+1):
    #     idx = int(numpy.floor((i / (num_ticks)) * (pole_freqs_size-1)))
    #     ticks.append(str(int(pole_freqs[idx])))
    #
    # plt.yticks(test, ticks)
    #
    # plt.ylim(numpy.min(spike_ids),numpy.max(spike_ids))
    # plt.xlim(0,numpy.max(scaled_times))
    # plt.xlabel('time (s)')
    # plt.ylabel('AN fibre best frequency (Hz)')

    #numpy.save("./spike_trains_asc_test_{}s.npy".format(int(duration)),scaled_times)
    numpy.savez_compressed('./spinnakear_asc_des_{}s_{}dB'.format(int(duration),int(dBSPL)), scaled_times=scaled_times,
             onset_times=onset_times,dBSPL=dBSPL)

#plot_output_spikes(spike_trains,plotter=plt,markersize=1,color='black')'''

else:
    plt.figure()
    drnl_time = numpy.linspace(0,len(audio_data)/Fs,len(audio_data))
   # plt.plot(drnl_time,drnl.T)
    plt.plot(drnl[1][:])
    #plt.ylim((-0.5e-7,0.5e-7))

    #plt.figure()
    #plt.plot(drnl[-1][:])
    #plt.ylim((-0.5e-7,0.5e-7))
  #  plt.xlim((0,len(audio_data)/Fs))
    print "max_lsr={}".format(numpy.max(drnl[0][:]))
    print "max_hsr={}".format(numpy.max(drnl[1][:]))

print "total stimulus time=",(len(audio_data)/Fs)
# plt.show()

#print "audio_data[9998]{:.100e}".format(audio_data[9998])
#print "single audio_data[9998]{:.100e}".format(numpy.float32(audio_data[9998]))

# Save the results
#numpy.save("./results.npy",drnl[1][:])
#savemat("/home/rjames/Dropbox (The University of Manchester)/EarProject/MAP_BS/MAP/results.mat",mdict={'spinn':drnl[1][:]},long_field_names=True)
#numpy.save("/home/rjames/Dropbox (The University of Manchester)/EarProject/spike_trains_no5.npy", [spike_trains,duration,Fs])
#numpy.save("/home/rjames/Dropbox (The University of Manchester)/EarProject/spike_trains_10sp_2num_5rep.npy", [spike_trains,duration,Fs])
#numpy.savetxt("./results.csv",drnl, fmt="%e", delimiter=",")

#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/results.csv", drnl[1][:], fmt="%.1000e", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/results.csv", audio_data, fmt="%.38e", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/complete.txt",[1],fmt="%f")
