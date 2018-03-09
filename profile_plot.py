import pylab as plt
import numpy as np
from profile_average import profile_average
from glob import glob



#profile_ome = np.load('./0_0_1_profile.npy')
#profile_drnl = np.load('./0_0_2_profile.npy')
#profile_ihc = np.load('./0_0_9_profile.npy')

def profile_plot():
    #profile_average(2048)
    master_profile = np.load('./master_profiles.npy')
    num_segs = 13152./96.
    channels =[2**(i+1) for i in range(11)]
    channels.append(3000.)
    omes = master_profile[1::4]
    omes = [i * num_segs * 1e-3 for i in omes]
    drnls = master_profile[2::4]
    drnls = [i * num_segs * 1e-3 for i in drnls]
    ihcs = master_profile[3::4]
    ihcs = [i * num_segs * 1e-3 for i in ihcs]

    rt = 13152./22050.
    #plt.figure()
    #plt.plot(channels,omes)
    #plt.plot(channels,drnls)
    #plt.plot(channels,ihcs)

    complete = np.asarray([omes,drnls,ihcs])
    complete = complete.T

    np.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/SpiNNakEar_profiles_3.csv", complete, fmt="%e", delimiter=",")

def single_profile_plot():

    profile_drnls=[]
    profile_ihcans=[]
    fnames = glob('./*_profile.npy')

    for f in fnames:
        placement=f.lstrip('./')
        placement=placement.rstrip('_profile.npy')
        if placement[0] == '0' and placement[2]=='0':#0 chip contains ome
            if placement[-1] == '1' and placement[-2]=='_':#ome
                profile_ome = np.load(f)
            elif placement[4] == '2' or placement[4] == '8':
                profile_drnls.append(np.load(f))
            else:
                profile_ihcans.append(np.load(f))
        elif (placement[-1] == '1' and placement[-2] == '_') or placement[-1] == '7':#any other chip drnl location shifts
            profile_drnls.append(np.load(f))
        else:#ihcan
            profile_ihcans.append(np.load(f))

    #profile_ihc = np.load('./0_0_9_profile.npy')
    #ihcs = [i * num_segs * 1e-3 for i in profile_ihc]

    #plt.plot(profile_ome)
    min_len = len(min(profile_drnls,key=len))
    profile_drnls = [d[:min_len] for d in profile_drnls]
    #plt.plot(np.array(profile_drnls).T)
    min_len = len(min(profile_ihcans,key=len))
    profile_ihcans = [i[:min_len] for i in profile_ihcans]
    plt.plot(np.array(profile_ihcans).T)

    #plt.plot(profile_drnl+profile_ome)
    #plt.plot(profile_ihc+profile_drnl+profile_ome)

    #print "max ihc execution time: {}ms".format(np.max(profile_ihc))

    #print "max total execution time: {}ms".format(np.max(profile_ihc+profile_drnl+profile_ome))"""

    plt.show()

single_profile_plot()