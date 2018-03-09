import pylab as plt
import numpy as np
from glob import glob
#find 0 chip (contains OME)


def profile_average(n_ch,plt=None):
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

    profile_ihcans = [p for p in profile_ihcans if len(p)==137]
    profile_drnls = [p for p in profile_drnls if len(p)==137]
    profile_drnl = np.mean(profile_drnls,0)
    profile_ihc = np.mean(profile_ihcans,0)

    if plt!=None:
        plt.plot(profile_ome)
        plt.plot(profile_drnl)
        #plt.plot(profile_drnl+profile_ome)
        plt.plot(profile_ihc)
        #plt.plot(profile_ihc+profile_drnl+profile_ome)
        plt.xlim((0,len(profile_ome)))
        #print "max ihc execution time: {}ms".format(np.max(profile_ihc))

        #print "max total execution time: {}ms".format(np.max(profile_ihc+profile_drnl+profile_ome))
        plt.show()


    profiles=np.asarray([n_ch,np.mean(profile_ome),np.mean(profile_drnl),np.mean(profile_ihc)])
    master_profiles = np.load("./master_profiles.npy")
    master_profiles = master_profiles[0:44]
    master_profiles = np.append(master_profiles,profiles)
    np.save("./master_profiles.npy",master_profiles)


    #np.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/profile{}.csv".format(n_ch),
    #           profiles, fmt="%e", delimiter=",")

