import pylab as plt
import numpy as np

profile_ome = np.load('./0_0_1_profile.npy')
profile_drnl = np.load('./0_0_2_profile.npy')
profile_ihc = np.load('./0_0_9_profile.npy')

plt.plot(profile_ome)
plt.plot(profile_drnl)
#plt.plot(profile_drnl+profile_ome)
plt.plot(profile_ihc)
#plt.plot(profile_ihc+profile_drnl+profile_ome)

print "max ihc execution time: {}ms".format(np.max(profile_ihc))

#print "max total execution time: {}ms".format(np.max(profile_ihc+profile_drnl+profile_ome))
plt.show()