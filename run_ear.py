import numpy
import model_launch_framework
#from mcmc_examples.lighthouse.lighthouse_model import LightHouseModel
from scipy.io import wavfile
import pylab as plt


#fs, audio_data=wavfile.read('4kHz_40dB.wav')

audio_data=numpy.fromfile("./c_models/load_files/load1_1",dtype='float32')
#audio_data=numpy.fromfile("./c_models/load_files/load1_1kate22k",dtype='float32')
pole_freqs=numpy.fromfile("./c_models/load_files/pole_freqs",dtype='float32')
#pole_freqs=[500, 1000, 2000, 3000, 4000, 5000, 6000, 7000]
pole_freqs=[4000,4000]
#plt.figure()
#plt.plot(audio_data[12000:15000])#[12000:20500])
plt.figure()
plt.plot(audio_data)
#plt.show()


#TODO: create auditory model vertex instances

#create framework of connected model vertices and run
samples = model_launch_framework.run_model(
    audio_data, n_chips=len(pole_freqs)/2,n_drnl=2,pole_freqs=pole_freqs,n_ihcan=5,fs=22050,resample_factor=1)

# Save the results
numpy.save("results.npy", samples)
numpy.savetxt("results.csv", samples, fmt="%f", delimiter=",")


