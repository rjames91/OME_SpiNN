import numpy
import model_launch_framework_binaural
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
tone_13 = generate_signal(freq=13500,dBSPL=dBSPL,duration=1.,
                       modulation_freq=0.,fs=Fs,ramp_duration=0.0025,plt=None,silence=True)
tone_1 = generate_signal(freq=1000,dBSPL=dBSPL,duration=1.,
                       modulation_freq=0.,fs=Fs,ramp_duration=0.0025,plt=None,silence=True)
matches_l = generate_signal(signal_type='file',dBSPL=dBSPL,fs=Fs,ramp_duration=0.0025,silence=True,
                            file_name='./binaural_matches_6s.wav',plt=None,channel=0)
matches_r = generate_signal(signal_type='file',dBSPL=dBSPL,fs=Fs,ramp_duration=0.0025,silence=True,
                            file_name='./binaural_matches_6s.wav',plt=None,channel=1)
matches = numpy.asarray([matches_l,matches_r])

sounds = [matches,tone_13,tone_1]#[matches_l]#[asc,des,a,i,u]#
# sounds = [numpy.zeros((2,int(Fs)))]#[matches_l]#[asc,des,a,i,u]#
# check if any sounds are in stereo
num_channels=1
for sound in sounds:
    if len(sound.shape)>1:
        num_channels=2
audio_data = [[] for _ in range(num_channels)]
required_total_time = 60.
onset_times = [[[]for _ in range(num_channels)]for _ in range(len(sounds))]

while 1:
    rand_choice = numpy.random.randint(len(sounds))
    chosen_samples = sounds[rand_choice]
    if len(chosen_samples.shape) > 1:#stereo
        for i,channel in enumerate(chosen_samples):
            onset_found = False
            for sample in channel:
                if not onset_found and abs(sample) > 1e-10:
                    onset_times[rand_choice][i].append(numpy.round(1000. * (len(audio_data[i]) / Fs)))
                    onset_found = True
                audio_data[i].append(sample)

    else:
        onset_found = False
        for sample in chosen_samples:
            if not onset_found and abs(sample)>1e-10:
                onset_times[rand_choice][0].append(numpy.round(1000.*(len(audio_data[0])/Fs)))
                onset_found = True
            audio_data[0].append(sample)
            #add zero input to other ear
            if num_channels > 1:
                audio_data[1].append(0.)

    if int(len(audio_data[0]) / Fs) > required_total_time:
        break

if Fs > 34000.:
    pole_freqs = numpy.logspace(1.477,4.25,100)#full hearing spectrum
    rt = False
else:
    pole_freqs = numpy.logspace(1.477,3.95,100)#up to 9k range
    rt= True

binaural_audio_data = []
#check audio data can be divided evenly into 100 sample segements
for ear in range(num_channels):
    binaural_audio_data.append(numpy.asarray(audio_data[ear][0:int(numpy.floor(len(audio_data[ear])/seg_size)*seg_size)]))

binaural_audio_data = numpy.asarray(binaural_audio_data)

duration = binaural_audio_data[0].size/Fs

#create framework of connected model vertices and run
samples = model_launch_framework_binaural.run_model(
    binaural_audio_data, n_chips=numpy.ceil(len(pole_freqs)/2)*num_channels,n_drnl=2,
    pole_freqs=pole_freqs,n_ihcan=5,fs=Fs,resample_factor=resample_factor,num_macks=4,
    bitfield=bitfield,rt=rt,profile=profile)

binaural_scaled_times = []
binaural_drnl = []
#begin for each ear

for ear in range(num_channels):
    #convert to spike train
    spike_trains=[]
    spike_times = []
    for _ in range(len(pole_freqs)*10):#TODO add n_fibres per IHC parameter to replace 10
        spike_times.append([])
    spike_index=0
    ihc_index=0

    #obtain list of IHCAN outputs
    ihc_output = [samples[ear][x:x+2*int(numpy.floor(len(binaural_audio_data[ear])/seg_size))*seg_size]
                  for x in xrange(0, len(samples[ear]), 2*int(numpy.floor(len(binaural_audio_data[ear])/seg_size))*seg_size)]
    drnl=numpy.zeros((len(pole_freqs)*10,int(numpy.floor(len(binaural_audio_data[ear])/seg_size))*seg_size),dtype=numpy.float32)
    for ihc in ihc_output:
        #obtain fibre response
        hsr = [ihc[x:x+seg_size] for x in xrange(0,len(ihc),seg_size*2)]
        HSR=numpy.concatenate(hsr,axis=0)
        drnl[-spike_index][:]=HSR
        if bitfield:
            #find non zero indicies for fibre and record in spike trains list
            idxs = numpy.nonzero(HSR)
            spike_times[spike_index]=idxs[0]

        #increment spike_index
        spike_index += 1

        #obtain fibre response
        lsr=[ihc[x:x+seg_size] for x in xrange(seg_size,len(ihc),seg_size*2)]
        LSR = numpy.concatenate(lsr,axis=0)
        drnl[-spike_index][:]=LSR
        if bitfield:
            # find non zero indicies for fibre and record in spike trains list
            idxs = numpy.nonzero(LSR)
            spike_times[spike_index]=idxs[0]

        # increment spike_index
        spike_index += 1
        # increment ihc_index
        ihc_index += 1

    if bitfield:
        max_time = 0
        for neuron in spike_times:
            if neuron.size > 0:
                if neuron.max()>max_time:
                    max_time = neuron.max()
        scale_factor = duration/max_time
        scaled_times = [1000 * spike_time * scale_factor for spike_time in spike_times]
        spike_raster_plot_8(scaled_times, plt, duration, len(pole_freqs)*10 + 1, 0.001, title="pre pop activity ear {}".format(ear))
        binaural_scaled_times.append(scaled_times)

    else:
        plt.figure()
        drnl_time = numpy.linspace(0,len(binaural_audio_data[ear])/Fs,len(binaural_audio_data[ear]))
        plt.plot(drnl[1][:])
        print "max_lsr={}".format(numpy.max(drnl[0][:]))
        print "max_hsr={}".format(numpy.max(drnl[1][:]))
        binaural_drnl.append(drnl)

print "total stimulus time = {}".format(duration)

# Save the results
if bitfield:
    numpy.savez_compressed('./spinnakear_matches_{}s_{}dB'.format(int(duration),int(dBSPL)), scaled_times=binaural_scaled_times,
             onset_times=onset_times,dBSPL=dBSPL)
else:
    numpy.savez_compressed('./spinnakear_13.5_1_kHz_{}s_{}dB_vrr'.format(int(duration), int(dBSPL)),
                           drnl=binaural_drnl,onset_times=onset_times, dBSPL=dBSPL)
#numpy.save("./results.npy",drnl[1][:])
#savemat("/home/rjames/Dropbox (The University of Manchester)/EarProject/MAP_BS/MAP/results.mat",mdict={'spinn':drnl[1][:]},long_field_names=True)
#numpy.save("/home/rjames/Dropbox (The University of Manchester)/EarProject/spike_trains_no5.npy", [spike_trains,duration,Fs])
#numpy.save("/home/rjames/Dropbox (The University of Manchester)/EarProject/spike_trains_10sp_2num_5rep.npy", [spike_trains,duration,Fs])
#numpy.savetxt("./results.csv",drnl, fmt="%e", delimiter=",")

#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/results.csv", drnl[1][:], fmt="%.1000e", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/results.csv", audio_data, fmt="%.38e", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/complete.txt",[1],fmt="%f")
plt.show()