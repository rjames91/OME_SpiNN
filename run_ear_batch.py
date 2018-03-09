import numpy
import model_launch_framework
from scipy.io import wavfile
import pylab as plt
import math
from vision.spike_tools.vis.vis_tools import plot_output_spikes
from signal_prep import *
from profile_average import profile_average
from profile_plot import profile_plot
import os

master_profiles = []
numpy.save("./master_profiles.npy", master_profiles)

tests=[2**(i+1) for i in range(11)]
tests.append(3000)
#tests = [2048,3000]

for channels in tests:
    os.system("rm ./*_profile.npy")
    Fs = 22050.#24000.#44100.#40000.#100000.#
    seg_size = 96
    bitfield = True
    #audio_data=numpy.fromfile("./c_models/load_files/load1_1",dtype='float32')
    #audio_data=numpy.fromfile("./c_models/load_files/load1_1_6k_22k",dtype='float32')
    #audio_data=numpy.fromfile("./c_models/load_files/load1_1_6k_44k",dtype='float32')
    #audio_data=numpy.fromfile("./c_models/load_files/load1_1_chirp",dtype='float32')
    #audio_data=numpy.fromfile("./c_models/load_files/load1_1kate_22k",dtype='float32')
    #audio_data=numpy.fromfile("./c_models/load_files/load1_1vowels_22k",dtype='float32')
    #audio_data = wavfile.read('./vowels_44_1kHz.wav')
    #audio_data = numpy.float32(audio_data[1])
    #audio_data = generate_signal(freq=6900,dBSPL=68.,duration=0.4,
    #                             modulation_freq=0.,fs=Fs,ramp_duration=0.01,plt=plt,silence=False)

    audio_data = generate_signal(signal_type="sweep_tone",freq=[30,8000],dBSPL=50.,duration=0.5,
                                 modulation_freq=0.,fs=Fs,ramp_duration=0.0025,plt=None,silence=True)
    #audio_data = numpy.concatenate((audio_data,audio_data))

    #audio_data = generate_signal(signal_type='file',dBSPL=60.,fs=Fs,ramp_duration=0.01,
    #                             file_name='./vowels_a.wav',plt=plt)

    #numpy.save('../Brainstem/audio_data.npy',audio_data)
    #concha = test_filter(audio_data,0.1783,0,-0.1783,1,-1.3477,0.6433)
    #numpy.savetxt("./concha.csv",concha, fmt="%e", delimiter=",")
    #audio_data=numpy.fromfile("./c_models/load_files/load1_1",dtype='float32')
    #pole_freqs=numpy.fromfile("./c_models/load_files/pole_freqs_125",dtype='float32')

    pole_freqs = numpy.logspace(2,3.95,channels)

    #pole_freqs=numpy.empty(71)#TODO:investigate intermitent simulation failures
    #pole_freqs.fill(6900.0)
    #pole_freqs=numpy.empty(16)
    #pole_freqs.fill(4000.0)
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
    audio_data = audio_data[0:int(numpy.floor(len(audio_data)/seg_size)*seg_size)]
    #plt.figure()
    #plt.plot(audio_data)
    #plt.show

    #create framework of connected model vertices and run
    samples = model_launch_framework.run_model(
        audio_data, n_chips=numpy.ceil(len(pole_freqs)/2),n_drnl=2,
        pole_freqs=pole_freqs,n_ihcan=5,fs=Fs,resample_factor=1,num_macks=4,
        bitfield=bitfield,rt=False)

    profile_average(channels,plt=None)

#numpy.save("./samples.npy",samples)

#samples=numpy.load("./samples.npy")

#call profile plot
profile_plot()

# Save the results
#numpy.save("/home/rjames/Dropbox (The University of Manchester)/EarProject/spike_trains_6k_640fib_50dB.npy", spike_trains)
#numpy.savetxt("./results.csv",drnl, fmt="%e", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/results.csv", samples, fmt="%e", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/complete.txt",[1],fmt="%f")
