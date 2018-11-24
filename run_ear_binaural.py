import numpy
import model_launch_framework_binaural
from scipy.io import wavfile
import pylab as plt
import math
from signal_prep import *
from scipy.io import savemat, loadmat

Fs = 44100.#22050.#24000.#34000.#10000-0.#40000.#
seg_size = 96
bitfield = True#False#
profile = False#True#
psth = False
resample_factor = 1.
dBSPL = 30.
sweep_duration = 0.1

asc = generate_signal(signal_type="sweep_tone", freq=[30, 8000], dBSPL=dBSPL, duration=sweep_duration,
                            modulation_freq=0., fs=Fs, ramp_duration=0.0025, plt=None, silence=True,ascending=True)
des = generate_signal(signal_type="sweep_tone", freq=[30, 8000], dBSPL=dBSPL, duration=sweep_duration,
                            modulation_freq=0., fs=Fs, ramp_duration=0.0025, plt=None, silence=True,ascending=False)
yes = generate_signal(signal_type='file',dBSPL=dBSPL,fs=Fs,ramp_duration=0.0025,silence=True,
                            file_name='./yes_samples/yes_edit1.wav',plt=None)
a = generate_signal(signal_type='file',dBSPL=dBSPL,fs=Fs,ramp_duration=0.0025,silence=True,
                            file_name='./a.wav',plt=None)
i = generate_signal(signal_type='file',dBSPL=dBSPL,fs=Fs,ramp_duration=0.0025,silence=True,
                            file_name='./i.wav',plt=None)
u = generate_signal(signal_type='file',dBSPL=dBSPL,fs=Fs,ramp_duration=0.0025,silence=True,
                            file_name='./u.wav',plt=None)
tone_13 = generate_signal(freq=13500,dBSPL=dBSPL,duration=0.05,
                       modulation_freq=0.,fs=Fs,ramp_duration=0.0025,plt=None,silence=True)
tone_1 = generate_signal(freq=1000,dBSPL=dBSPL,duration=0.05,
                       modulation_freq=0.,fs=Fs,ramp_duration=0.0025,plt=None,silence=True)
tone_1_quiet = generate_signal(freq=1000,dBSPL=dBSPL/2.,duration=0.05,
                       modulation_freq=0.,fs=Fs,ramp_duration=0.0025,plt=None,silence=True)
matches_l = generate_signal(signal_type='file',dBSPL=dBSPL,fs=Fs,ramp_duration=0.0025,silence=True,
                            file_name='./binaural_matches_6s.wav',plt=None,channel=0)
matches_r = generate_signal(signal_type='file',dBSPL=dBSPL,fs=Fs,ramp_duration=0.0025,silence=True,
                            file_name='./binaural_matches_6s.wav',plt=None,channel=1)
matches = numpy.asarray([matches_l,matches_r])

stereo_1k = numpy.asarray([tone_1,tone_1_quiet])

sounds_dict = { "matches":matches,
                "matches_l":matches_l,
                "matches_r":matches_r,
                "tone_1":tone_1,
                "stereo_1k":stereo_1k,
                "tone_13":tone_13,
                "a":a,
                "i":i,
                "u":u,
                "asc":asc,
                "des":des,
                "yes":yes
}

#choose test stimuli here
stimulus_list = ['tone_1']
# check if any stimuli are in stereo
num_channels=1
for sound_string in stimulus_list:
    sound = sounds_dict[sound_string]
    if len(sound.shape)>1:
        num_channels=2
audio_data = [[] for _ in range(num_channels)]
# for channel in audio_data:
#     for _ in range(100):
#         channel.append(0.)

required_total_time = 1.
onset_times = [[[]for _ in range(num_channels)]for _ in range(len(stimulus_list))]

chosen_stimulus_list=[]
while 1:
    rand_choice = numpy.random.randint(len(stimulus_list))
    chosen_string = stimulus_list[rand_choice]
    if chosen_string not in chosen_stimulus_list:
        chosen_stimulus_list.append(chosen_string)
    chosen_samples = sounds_dict[chosen_string]
    if len(chosen_samples.shape) > 1:#stereo
        for i,channel in enumerate(chosen_samples):
            onset_found = False
            for sample in channel:
                if not onset_found and abs(sample) > 1e-10:
                    #onset time, duration tuple (both in ms)
                    onset_times[rand_choice][i].append(numpy.round(1000. * (len(audio_data[i]) / Fs)))
                    onset_found = True
                audio_data[i].append(sample)

    else:
        onset_found = False
        for sample in chosen_samples:
            if not onset_found and abs(sample) > 1e-10:
                onset_times[rand_choice][0].append(numpy.round(1000.*(len(audio_data[0])/Fs)))
                onset_found = True
            audio_data[0].append(sample)
            #add zero input to other ear
        if num_channels > 1:
            # silence_amp = 1. * 28e-6 * 10. ** ((dBSPL-60.) / 20.)
            silence_amp = 1. * 28e-6 * 10. ** (-20. / 20.) # -20dBSPL silence noise
            silence_samples = numpy.random.rand(chosen_samples.size) * silence_amp
            for sample in silence_samples:
                audio_data[1].append(sample)

    if len(audio_data[0]) / Fs > required_total_time:
        break

if Fs > 34000.:
    pole_freqs = numpy.logspace(1.477,4.25,100)#full hearing spectrum
    rt = False
else:
    pole_freqs = numpy.logspace(1.477,3.95,100)#up to 9k range
    rt = True

binaural_audio_data = []
#check audio data can be divided evenly into 100 sample segements
for ear in range(num_channels):
    binaural_audio_data.append(numpy.asarray(audio_data[ear][0:int(numpy.floor(len(audio_data[ear])/seg_size)*seg_size)]))

binaural_audio_data = numpy.asarray(binaural_audio_data)

duration = binaural_audio_data[0].size/Fs

#create framework of connected model vertices and run
samples,profiles = model_launch_framework_binaural.run_model(
    binaural_audio_data, n_chips=numpy.ceil(len(pole_freqs)/2)*num_channels,n_drnl=2,
    pole_freqs=pole_freqs,n_ihcan=5,fs=Fs,resample_factor=resample_factor,num_macks=4,
    bitfield=bitfield,rt=rt,profile=profile)

binaural_scaled_times = []
binaural_drnl = []

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
        if psth:
            psth_plot_8(plt, numpy.arange(1,len(scaled_times),2), scaled_times, bin_width=1. / 1000., duration=duration,
                        title="PSTH_AN_HSR ear {}".format(ear))
            psth_plot_8(plt, numpy.arange(0, len(scaled_times), 2), scaled_times, bin_width=1. / 1000., duration=duration,
                        title="PSTH_AN_LSR ear {}".format(ear))
        binaural_scaled_times.append(scaled_times)

    else:
        plt.figure("VRR output")
        drnl_time = numpy.linspace(0,len(binaural_audio_data[ear])/Fs,len(binaural_audio_data[ear]))
        x = numpy.linspace(0,1,len(drnl[1]))
        plt.plot(x,drnl[1][:])
        plt.plot(x,drnl[0][:])
        plt.xlim((0,1))
        print "max_lsr={}".format(numpy.max(drnl[0][:]))
        print "max_hsr={}".format(numpy.max(drnl[1][:]))
        binaural_drnl.append(drnl)

print "total stimulus time = {}".format(duration)

n_fibres = len(pole_freqs) * 10

# Save the results
results_directory = '/home/rjames/Dropbox (The University of Manchester)/EarProject/Pattern_recognition/spike_trains/IC_spikes'
stimulus_string = ""
for stim_string in chosen_stimulus_list:
    stimulus_string += stim_string + "_"
if bitfield:
    #TODO: add durations for each stimulus
    numpy.savez_compressed(results_directory+'/spinnakear_'+stimulus_string+'{}s_{}dB_{}fibres'.format(int(duration),int(dBSPL),int(n_fibres)),
                           scaled_times=binaural_scaled_times,onset_times=onset_times,dBSPL=dBSPL,profiles=profiles,Fs=Fs,audio_data=binaural_audio_data)
else:
    numpy.savez_compressed(results_directory+'/spinnakear_'+stimulus_string+'{}s_{}dB_non_spiking_{}fibres'.format(int(duration),int(dBSPL),int(n_fibres)),
                           drnl=binaural_drnl,onset_times=onset_times, dBSPL=dBSPL,Fs=Fs)

#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/results.csv", drnl[1][:], fmt="%.1000e", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/results.csv", audio_data, fmt="%.38e", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/complete.txt",[1],fmt="%f")

if profile:
    title_list = ["OME (ear{})","DRNL (ear{})","IHC (ear{})"]
    for ear in range(num_channels):
        for module in range(3):
            plt.figure(title_list[module].format(ear))
            # make sure profile readings are the same length
            if len(profiles[ear,module].shape)==1 and module>0:
                min_length = numpy.inf
                for readings in profiles[ear,module]:
                    if len(readings)<min_length:
                        min_length = len(readings)
                x = numpy.linspace(0, 1,min_length)
                mod_array = [reading[:min_length] for reading in profiles[ear,module]]
                plt.plot(mod_array)
            else:
                x = numpy.linspace(0, 1,profiles[ear,module].size)
                # plt.plot(x,profiles[ear,module])
                plt.plot(profiles[ear,module])
            plt.xlabel("segment index")
            plt.ylabel("time taken (ms)")

plt.figure("input audio")
for channel in binaural_audio_data:
    t = numpy.linspace(0,duration,len(channel))
    plt.plot(channel)
    plt.xlabel("time (s)")
plt.show()