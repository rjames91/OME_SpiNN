import numpy as np
from scipy.io import loadmat
import pylab as plt

spinn = np.load("./results.npy")

mat = loadmat("/home/rjames/Dropbox (The University of Manchester)/EarProject/MAP_BS/map_sv_64.mat")
sv_mat = mat["map_output_sv"][0]

mat  =loadmat("/home/rjames/Dropbox (The University of Manchester)/EarProject/MAP_BS/map_drnl.mat")
drnl_mat = mat["map_output"][0]

mat  =loadmat("/home/rjames/Dropbox (The University of Manchester)/EarProject/MAP_BS/map_vrr.mat")
vrr_mat = mat["map_output"][0]

mat = loadmat("/home/rjames/Dropbox (The University of Manchester)/EarProject/MAP_BS/map_coeffs.mat")
b = mat["pressure2dispb"][0]
a = mat["pressure2dispa"][0]

matlab = loadmat("/home/rjames/Dropbox (The University of Manchester)/EarProject/MAP_BS/map_64.mat")

MAP = matlab["map_output"][0]

matlab = loadmat("/home/rjames/Dropbox (The University of Manchester)/EarProject/MAP_BS/map.mat")
MAP_single = matlab["map_output"][0]

spinn = spinn[96:]
#MAP_32 =MAP[:len(spinn)].astype(np.float32)#sv_mat[:len(spinn)].astype(np.float32)#
MAP_32 = vrr_mat[:len(spinn)]#.astype(np.float32)#drnl_mat[:len(spinn)].astype(np.float32)#
diff = np.subtract(MAP_32[2:],spinn[2:])

#stapes displacement filter coeffs
a = -0.99937168146928201384326939660240896046161651611328125
b = 1.0 + a

spinn_diff_10181 = 4.9550076008674011441657644452083651044767796667311898772823042236268520355224609375000000000000000000e-26

if np.max(np.abs(diff))>0:
    first_non_zero = np.nonzero(diff)[0].min()
    print "first non-zero diff:{} (diff = {:.100e}".format(first_non_zero+2,diff[first_non_zero])
    print "Matlab[{}](as float) = {:.100e}".format(first_non_zero+2,MAP_32[first_non_zero+2])
    print "SpiNNaker[{}] = {:.100e}".format(first_non_zero+2,spinn[first_non_zero+2])
    python_check_eq = sv_mat[first_non_zero]*b - a*(MAP[first_non_zero-1]+spinn_diff_10181)
    # print "filter equation check(as float) = {:.100e}".format(np.float32(python_check_eq))
    # print "Matlab [{}](as double) = {:.100e}".format(first_non_zero-1,MAP[first_non_zero-1])
    print
    max_non_zero = np.argmax(np.abs(diff))
    print "max non-zero diff:{} (diff = {:.100e}".format(max_non_zero+2,diff[max_non_zero])
    print "Matlab[{}](as float) = {:.100e}".format(max_non_zero+2,MAP_32[max_non_zero+2])
    print "SpiNNaker[{}] = {:.100e}".format(max_non_zero+2,spinn[max_non_zero+2])
    python_check_eq = sv_mat[max_non_zero] * b - a * MAP[max_non_zero - 1]
    #print "filter equation check(as float) = {:.100e}".format(np.float32(python_check_eq))
#nonz = np.nonzero(diff)

#print "SpiNNaker[{}] = {:.100e}".format(10181,spinn[10181])

#print "{:.100e}".format(stapes_velocity[10196])
#print "{:.100e}".format(MAP[10194])

#test = spinn[999]#stapes_velocity[10195]*b-MAP[10194]*a#[1]#MAP[10194]#
#print "sample 9998 - -1.112e-8={:.100e}".format(spinn[9998])
#print "sample 9999 - -1.112e-8={:.100e}".format(spinn[9999])
# print "sample 9996={:.100e}".format(spinn[9996])
# print "sample 9997={:.100e}".format(spinn[9997])
# print "sample 9998={:.100e}".format(spinn[9998])
# print "sample 9999={:.100e}".format(spinn[9999])

print""
#diff_test = test - MAP[10196].astype(np.float32)
#print "diff={}".format(diff_test)
prop = np.abs(diff) / MAP_32[2:]
error = 100.*(np.nan_to_num(prop))

plt.figure()
plt.plot(diff)
plt.plot(error)
# plt.plot(MAP_32)

#plt.plot(spinn)
plt.title("difference between Matlab and SpiNNaker")
#plt.ylabel("stapes displacement")
#plt.xlabel("sample index")
plt.show()