import matplotlib.pylab as plt
import numpy as np
import matplotlib.ticker as ticker
from signal_prep import *

[spike_trains,scale_factor]=np.load('./spike_trains.npy')
spike_ids = [neuron_id for (neuron_id, spike_time) in spike_trains]

lsr = [(neuron_id,spike_time) for (neuron_id,spike_time) in spike_trains if neuron_id%2==0]
hsr = [(neuron_id,spike_time) for (neuron_id,spike_time) in spike_trains if neuron_id%2!=0]

spike_times = [spike_time for (neuron_id, spike_time) in lsr]
scaled_times = [spike_time * scale_factor for spike_time in spike_times]
spike_ids = [neuron_id for (neuron_id, spike_time) in lsr]
spike_ids[:] = [neuron_id + 1 for neuron_id in spike_ids]
max_id = np.max(spike_ids)

spike_times = [spike_time for (neuron_id, spike_time) in hsr]
h_scaled = [spike_time * scale_factor for spike_time in spike_times]
scaled_times+=h_scaled
s_ids = [neuron_id + max_id for (neuron_id, spike_time) in hsr]
spike_ids+=s_ids

##plot results
plt.figure()
plt.plot(scaled_times, spike_ids, '.', markersize=1,
         markerfacecolor='black', markeredgecolor='none',
         markeredgewidth=0)

ax=plt.gca()
pole_freqs = np.logspace(1.477,4.25,3000)
#pole_freqs =np.concatenate((pole_freqs,pole_freqs))
#max_id = np.max(spike_ids)
size = len(spike_ids)
num_ticks = 5.#16.0
increment = int(size/(num_ticks-1))
#test = spike_ids[::increment]
pole_freqs_size = len(pole_freqs)
test = [spike_ids[0],max_id-800,max_id+800,spike_ids[-1]]
#increment = int(pole_freqs_size/(num_ticks-1))-1
#ticks = pole_freqs[::increment]
ticks = [pole_freqs[0],pole_freqs[-1],pole_freqs[0],pole_freqs[-1]]
ticks_str= [str(int(np.round(i))) for i in ticks]

plt.yticks(test, ticks_str)

r_ticks = [max_id/2,max_id+max_id/2]
labels = ['LSR','HSR']
ax.yaxis.set_minor_locator(ticker.FixedLocator(r_ticks))
ax.yaxis.set_minor_formatter(ticker.FixedFormatter(labels))
ax.tick_params(axis="y", which="minor", direction="out",
                       left=0, right=1, labelleft=0, labelright=1)
"""
for tick in ax.yaxis.get_minorticklabels():
    if tick._text == 'HSR' or tick._text == 'LSR':
        tick.set_weight('bold')
        print
"""
plt.ylim(np.min(spike_ids), np.max(spike_ids))
plt.xlim(0, np.max(scaled_times))
plt.xlabel('time (s)')
plt.ylabel('AN fibre best frequency (Hz)')
plt.title('AN - SpiNNaker')

#PSTH_AN = generate_psth(numpy.arange(800,900),spike_trains,0.001,
#                     duration,scale_factor=an_scale_factor,Fs=Fs)
#x = numpy.arange(0,duration,duration/float(len(PSTH_AN)))
#plt.figure('PSTH_AN')
#plt.plot(x,PSTH_AN)

plt.show()