from signal_prep import *
import numpy as np
import matplotlib.pylab as plt
input_directory = '/home/rjames/Dropbox (The University of Manchester)/EarProject/Pattern_recognition/spike_trains/IC_spikes'

# ear_file = np.load(input_directory+'/spinnakear_1_kHz_2s_30dB.npz')
ear_file = np.load(input_directory+'/spinnakear_13.5_1_kHz_1s_-20dB_ome.npz')
ear_file_round = np.load(input_directory+'/spinnakear_13.5_1_kHz_1s_-20dB_ome_round.npz')
bin_drnl_round = ear_file_round['drnl']
drnl_round = bin_drnl_round[0]
bin_drnl = ear_file['drnl']
drnl=bin_drnl[0]
x = numpy.linspace(0, 1, len(drnl[1]))
plt.figure("OME output original")
plt.plot(x, drnl[1][:])
plt.plot(x, drnl_round[1][:])
plt.figure()
diff = (np.abs(drnl[1]-drnl_round[1]))
# error = 100./(np.abs(drnl[1]-drnl_round[1])/drnl[1])
plt.plot(x,diff)

# plt.xlim((0,len(drnl[1])))
plt.xlim((0, 1))
plt.show()

# profiles = ear_file['profiles']
#
# if profiles.shape[0]>1:
#     n_ears = 2
# else:
#     n_ears = 1
#
# title_list = ["OME (ear{})","DRNL (ear{})","IHC (ear{})"]
# for ear in range(n_ears):
#     for module in range(3):
#         plt.figure(title_list[module].format(ear))
#         # make sure profile readings are the same length
#         if len(profiles[ear,module].shape)==1 and module>0:
#             min_length = np.inf
#             for readings in profiles[ear,module]:
#                 if len(readings)<min_length:
#                     min_length = len(readings)
#             mod_array = [reading[:min_length] for reading in profiles[ear,module]]
#             plt.plot(mod_array)
#         else:
#             plt.plot(profiles[ear,module])
#         plt.xlabel("segment index")
#         plt.ylabel("time taken (ms)")
#
# plt.show()
# onset_times = ear_file['onset_times']

# for ear in onset_times:
#     a=1

drnl = ear_file['drnl']

plt.figure()
plt.plot(drnl)