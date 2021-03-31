#!/usr/bin/env python3

import numpy as np
import h5py

""" Averages over num_events trento events and saves the profile as a numpy array to be read by run-hydro.py."""

num_events = 50000
grid_size = 500
centrality_percent = 0.2
num_central = int(num_events*centrality_percent)

profile = np.zeros([grid_size,grid_size])
b = 0
npart = 0
mult = 0
e2 = 0
e3 = 0
e4 = 0
e5 = 0

bvals = np.zeros(num_events)

with h5py.File('events.hdf','r') as f:
    i = 0
    for dset in f.values():
        bvals[i] = dset.attrs['b']
        i = i+1
        
bvals = np.sort(bvals)
bvalsallowed = bvals[:num_central]

with h5py.File('events.hdf','r') as f:
    for dset in f.values():
        if dset.attrs['b'] in bvalsallowed:  
            profile += np.array(dset)
            b += dset.attrs['b']
        # npart += dset.attrs['npart']
        # mult += dset.attrs['mult']
        # e2 += dset.attrs['e2']
        # e3 += dset.attrs['e3']
        # e4 += dset.attrs['e4']
        # e5 += dset.attrs['e5']
        
profile = profile/num_central
b = b/num_central
# npart = npart/100000
# mult = mult/100000
# e2 = e2/100000
# e3 = e3/100000
# e4 = e4/100000
# e5 = e5/100000

np.save('profiles/profilex15n13p0grid500central20.npy', profile)