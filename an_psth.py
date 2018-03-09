import numpy as np
import matplotlib.pylab as plt
from signal_prep import *

Fs = 38000.#22050.#

[spike_trains,scale_factor]=np.load('./spike_trains.npy')

an_scale_factor = 1./Fs#duration/numpy.max(spike_times)
spike_times = [spike_time for (neuron_id, spike_time) in spike_trains]

scaled_times = [spike_time * an_scale_factor for (id,spike_time) in spike_trains]
duration = np.max(scaled_times)

PSTH_AN = generate_psth(numpy.arange(1,719,2),spike_trains,bin_width=0.001,
                     duration=duration,scale_factor=an_scale_factor,Fs=Fs)
x = numpy.arange(0,duration,duration/float(len(PSTH_AN)))
plt.figure('PSTH_AN')
plt.plot(x,PSTH_AN)

plt.show()