#!/usr/bin/env python3

import numpy as np
import h5py

""" Averages over 100,000 trento events and saves the profile as a numpy array to be read by run-hydro.py.
For trento normalization factor of 13, 100000 Pb Pb events are first run and saved in events13.hdf. The profile is saved as profile13.npy. """

profile = np.zeros([300,300])
b = 0
npart = 0
mult = 0
e2 = 0
e3 = 0
e4 = 0
e5 = 0

bvals = np.zeros(100000)

# with h5py.File('events13.hdf','r') as f:
#     i = 0
#     for dset in f.values():
#         bvals[i] = dset.attrs['b']
#         i = i+1
        
# bvals = np.sort(bvals)

with h5py.File('events.hdf','r') as f:
    for dset in f.values():
        profile += np.array(dset)
        # b += dset.attrs['b']
        # npart += dset.attrs['npart']
        # mult += dset.attrs['mult']
        # e2 += dset.attrs['e2']
        # e3 += dset.attrs['e3']
        # e4 += dset.attrs['e4']
        # e5 += dset.attrs['e5']
        
profile = profile/100000
# b = b/100000
# npart = npart/100000
# mult = mult/100000
# e2 = e2/100000
# e3 = e3/100000
# e4 = e4/100000
# e5 = e5/100000

np.save('profilex15n13p1grid300.npy', profile)