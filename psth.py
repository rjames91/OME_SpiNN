from signal_prep import *
import numpy as np
import pylab as plt


spike_trains=numpy.load("/home/rjames/Dropbox (The University of Manchester)/EarProject/spike_trains_6k_640fib_50dB.npy")

spike_trains=numpy.load("./spike_trains.npy")
target_ids =range(0,710,2)
duration = (0.5*21984.)/22050.
psth=generate_psth(target_ids,spike_trains[0],bin_width=0.001,duration=duration,
                   scale_factor=spike_trains[1])
print len(spike_trains)

plt.figure()

x = numpy.arange(0,duration,duration/float(len(psth)))

#plt.bar(psth,x)
plt.plot(psth)


plt.show()