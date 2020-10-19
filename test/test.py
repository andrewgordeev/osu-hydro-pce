import freestream
import h5py
import numpy as np
import math
import subprocess
import frzout

"""
Stripped down version of hic-eventgen's run-events, to study the effects of freestreaming + hydro 
code on data from one trento event.
"""

grid_step = 0.075
grid_n = math.ceil(2*15/grid_step)
grid_max = 0.5*grid_n*grid_step

def initial_conditions():
    with h5py.File('../../trento/test/test.hdf', 'r') as f:
        for dset in f.values():
            ic = np.array(dset)
            yield ic

def run_hydro_pce(fs, event_size, coarse=False, dt_ratio=0.25):
    dxy = grid_step * (coarse or 1)
    ls = math.ceil(event_size/dxy)
    n = 2*ls + 1
    
    for fmt, f, arglist in [
            ('ed', fs.energy_density, [()]),
            ('u{}', fs.flow_velocity, [(1,), (2,)]),
            ('pi{}{}', fs.shear_tensor, [(1,1),(1,2),(2,2)]),
            ]:
        for a in arglist:
            X= f(*a)
            
            if coarse:
                X = X[::coarse, ::coarse]
                
            diff = X.shape[0] -n
            start = int(abs(diff)/2)
            
            if diff > 0:
                s = slice(start, start + n)
                X = X[s,s]
            elif diff < 0:
                Xn = np.zeros((n,n))
                s = slice(start, start + X.shape[0])
                Xn[s,s] = X
                X = Xn
            
            X.tofile(fmt.format(*a) + '.dat')
            
    dt = dxy * dt_ratio

    subprocess.run(['osu-hydro-pce', 't0=0.5', 'dt={}'.format(dt), 'dxy={}'.format(dxy), 'nls={}'.format(ls)],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    
    surface = np.fromfile('surface.dat', dtype='f8').reshape(-1,16)
  #  return surface            
            
    return dict(
        zip(['x', 'sigma', 'v'], np.hsplit(surface, [3,6,8])),
        pi=dict(zip(['xx','xy','yy'],surface.T[11:14])),
        Pi=surface.T[15])


def run_hydro(fs, event_size, coarse=False, dt_ratio=0.25):
    dxy = grid_step * (coarse or 1)
    ls = math.ceil(event_size/dxy)
    n = 2*ls + 1
    
    for fmt, f, arglist in [
            ('ed', fs.energy_density, [()]),
            ('u{}', fs.flow_velocity, [(1,), (2,)]),
            ('pi{}{}', fs.shear_tensor, [(1,1),(1,2),(2,2)]),
            ]:
        for a in arglist:
            X= f(*a)
            
            if coarse:
                X = X[::coarse, ::coarse]
                
            diff = X.shape[0] -n
            start = int(abs(diff)/2)
            
            if diff > 0:
                s = slice(start, start + n)
                X = X[s,s]
            elif diff < 0:
                Xn = np.zeros((n,n))
                s = slice(start, start + X.shape[0])
                Xn[s,s] = X
                X = Xn
            
            X.tofile(fmt.format(*a) + '.dat')
            
    dt = dxy * dt_ratio

    subprocess.run(['osu-hydro', 't0=0.5', 'dt={}'.format(dt), 'dxy={}'.format(dxy), 'nls={}'.format(ls)],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    
    surface = np.fromfile('surface.dat', dtype='f8').reshape(-1,16)
    return surface            
            
    # return dict(
    #     zip(['x', 'sigma', 'v'], np.hsplit(surface, [3,6,8])),
    #     pi=dict(zip(['xx','xy','yy'],surface.T[11:14])),
    #     Pi=surface.T[15])

def run_single_event(ic, pce=False):
    fs = freestream.FreeStreamer(ic, grid_max, 0.5)
    surface = run_hydro(fs, event_size=27, coarse=3)
    resdict = dict(zip(['x', 'sigma', 'v'], np.hsplit(surface, [3,6,8])))
    rmax = math.sqrt((resdict['x'][:,1:3]**2).sum(axis=1).max())
    
    if pce:
        results = run_hydro_pce(fs, event_size=rmax)
    else:
        results = run_hydro(fs, event_size=rmax)
        
    return results

initial_conditions = initial_conditions()

species = [
    ('pion', 211),
    ('kaon', 321),
    ('proton', 2212),
    ('Lambda', 3122),
    ('Sigma0', 3212),
    ('Xi', 3312),
    ('Omega', 3334),
]

    
float_t = '<f8'
int_t = '<i8'
complex_t = '<c16'

finalresults = np.empty((), dtype=[
        ('initial_entropy', float_t),
        ('nsamples', int_t),
        ('dNch_deta', float_t),
        ('dET_deta', float_t),
        ('dN_dy', [(s, float_t) for (s, _) in species]),
        ('mean_pT', [(s, float_t) for (s, _) in species]),
        ('pT_fluct', [('N', int_t), ('sum_pT', float_t), ('sum_pTsq', float_t)]),
        ('flow', [('N', int_t), ('Qn', complex_t, 8)]),
    ])


finalresultspce = np.empty((), dtype=[
        ('initial_entropy', float_t),
        ('nsamples', int_t),
        ('dNch_deta', float_t),
        ('dET_deta', float_t),
        ('dN_dy', [(s, float_t) for (s, _) in species]),
        ('mean_pT', [(s, float_t) for (s, _) in species]),
        ('pT_fluct', [('N', int_t), ('sum_pT', float_t), ('sum_pTsq', float_t)]),
        ('flow', [('N', int_t), ('Qn', complex_t, 8)]),
    ])


for n, ic in enumerate(initial_conditions, start=1):
    if n == 1:
        results = run_single_event(ic)
        resultsdict = dict(
        zip(['x', 'sigma', 'v'], np.hsplit(results, [3,6,8])),
        pi=dict(zip(['xx','xy','yy'],results.T[11:14])),
        Pi=results.T[15])
        resultspce = run_single_event(ic, pce=True)
        finalsurface = frzout.Surface(**resultsdict, ymax=2)
        finalsurfacepce = frzout.Surface(**resultspce, ymax = 2)
        finalresults['initial_entropy'] = ic.sum() * grid_step**2
        finalresultspce['initial_entropy'] = ic.sum() * grid_step**2
    else:
        continue
    

minsamples, maxsamples = 10, 1000  # reasonable range for nsamples
minparts = 10**5  # min number of particles to sample
nparts = 0  # for tracking total number of sampled particles
npartspce = 0

hrg_kwargs = dict(species='urqmd', res_width=True)
hrg = frzout.HRG(0.150, **hrg_kwargs)

for nsamples in range(1, maxsamples + 1):
    parts = frzout.sample(finalsurface, hrg)
    if parts.size == 0:
        continue
    nparts += parts.size
 #   print('#', parts.size, file=f)
    for p in parts:
        continue
 #       print(p['ID'], *p['x'], *p['p'], file=f)
    if nparts >= minparts and nsamples >= minsamples:
        break
    

for nsamplespce in range(1, maxsamples + 1):
    partspce = frzout.sample(finalsurfacepce, hrg)
    if partspce.size == 0:
        continue
    npartspce += partspce.size
 #   print('#', parts.size, file=f)
    for p in partspce:
        continue
 #       print(p['ID'], *p['x'], *p['p'], file=f)
    if npartspce >= minparts and nsamplespce >= minsamples:
        break    
    
finalresults['nsamples'] = nsamples
finalresultspce['nsamples'] = nsamplespce

charged = (parts['charge'] != 0)
abs_eta = np.fabs(parts['eta'])

finalresults['dNch_deta'] = \
    np.count_nonzero(charged & (abs_eta < .5)) / nsamples

ET_eta = .6
finalresults['dET_deta'] = \
    parts['ET'][abs_eta < ET_eta].sum() / (2*ET_eta) / nsamples

abs_ID = np.abs(parts['ID'])
midrapidity = (np.fabs(parts['y']) < .5)

pT = parts['pT']
phi = parts['phi']

for name, i in species:
    cut = (abs_ID == i) & midrapidity
    N = np.count_nonzero(cut)
    finalresults['dN_dy'][name] = N / nsamples
    finalresults['mean_pT'][name] = (0. if N == 0 else pT[cut].mean())

pT_alice = pT[charged & (abs_eta < .8) & (.15 < pT) & (pT < 2.)]
finalresults['pT_fluct']['N'] = pT_alice.size
finalresults['pT_fluct']['sum_pT'] = pT_alice.sum()
finalresults['pT_fluct']['sum_pTsq'] = np.inner(pT_alice, pT_alice)

phi_alice = phi[charged & (abs_eta < .8) & (.2 < pT) & (pT < 5.)]
finalresults['flow']['N'] = phi_alice.size
finalresults['flow']['Qn'] = [
    np.exp(1j*n*phi_alice).sum()
    for n in range(1, finalresults.dtype['flow']['Qn'].shape[0] + 1)
]

charged = (partspce['charge'] != 0)
abs_eta = np.fabs(partspce['eta'])

finalresultspce['dNch_deta'] = \
    np.count_nonzero(charged & (abs_eta < .5)) / nsamplespce

ET_eta = .6
finalresultspce['dET_deta'] = \
    partspce['ET'][abs_eta < ET_eta].sum() / (2*ET_eta) / nsamplespce

abs_ID = np.abs(partspce['ID'])
midrapidity = (np.fabs(partspce['y']) < .5)

pT = partspce['pT']
phi = partspce['phi']

for name, i in species:
    cut = (abs_ID == i) & midrapidity
    N = np.count_nonzero(cut)
    finalresultspce['dN_dy'][name] = N / nsamplespce
    finalresultspce['mean_pT'][name] = (0. if N == 0 else pT[cut].mean())

pT_alice = pT[charged & (abs_eta < .8) & (.15 < pT) & (pT < 2.)]
finalresultspce['pT_fluct']['N'] = pT_alice.size
finalresultspce['pT_fluct']['sum_pT'] = pT_alice.sum()
finalresultspce['pT_fluct']['sum_pTsq'] = np.inner(pT_alice, pT_alice)

phi_alice = phi[charged & (abs_eta < .8) & (.2 < pT) & (pT < 5.)]
finalresultspce['flow']['N'] = phi_alice.size
finalresultspce['flow']['Qn'] = [
    np.exp(1j*n*phi_alice).sum()
    for n in range(1, finalresultspce.dtype['flow']['Qn'].shape[0] + 1)
]

results_file = 'event.dat'
resultspce_file = 'eventpce.dat'
results_file.write(finalresults.tobytes())
resultspce_file.write(finalresultspce.tobytes())