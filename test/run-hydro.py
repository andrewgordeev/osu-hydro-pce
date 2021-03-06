import numpy as np
import math
import subprocess
import frzout
import freestream
try:
    from scipy.interpolate import CubicSpline
except ImportError:
    from scipy.interpolate import InterpolatedUnivariateSpline as CubicSpline
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import griddata

import matplotlib.pyplot as plt
    

"""
Generating plots 4-6 from arXiv:160406346v3 (Vovchenko electromagnetic probes paper). 
Much of the code is adapted from hic-eventgen's run-events at https://github.com/Duke-QCD/hic-eventgen
"""

grid_step = 30./300.
grid_n = math.ceil(2*15/grid_step)
grid_max = 0.5*grid_n*grid_step

def run_hydro(fs, event_size, coarse=False, dt_ratio=0.25, profile=None):
    """
    Run the initial condition contained in FreeStreamer object `fs` through
    osu-hydro-pce on a grid with approximate physical size `event_size` [fm].
    Return a dict of freeze-out surface data suitable for passing directly
    to frzout.Surface.

    Initial condition arrays are cropped or padded as necessary.

    If `coarse` is an integer > 1, use only every `coarse`th cell from the
    initial condition arrays (thus increasing the physical grid step size
    by a factor of `coarse`).  Ignore the user input `hydro_args` and
    instead run ideal hydro down to a low temperature.

    `dt_ratio` sets the timestep as a fraction of the spatial step
    (dt = dt_ratio * dxy).  The SHASTA algorithm requires dt_ratio < 1/2.
    """
    
    dxy = grid_step * (coarse or 1)
    ls = math.ceil(event_size/dxy)
    n = 2*ls + 1

    # for fmt, f, arglist in [
    #         ('ed', fs.energy_density, [()]),
    #         ('u{}', fs.flow_velocity, [(1,), (2,)]),
    #         ('pi{}{}', fs.shear_tensor, [(1, 1), (1, 2), (2, 2)]),
    # ]:
    #     for a in arglist:
    #         X = f(*a)

    #         if coarse:
    #             X = X[::coarse, ::coarse]

    #         diff = X.shape[0] - n
    #         start = int(abs(diff)/2)

    #         if diff > 0:
    #             # original grid is larger -> cut out middle square
    #             s = slice(start, start + n)
    #             X = X[s, s]
    #         elif diff < 0:
    #             # original grid is smaller
    #             #  -> create new array and place original grid in middle
    #             Xn = np.zeros((n, n))
    #             s = slice(start, start + X.shape[0])
    #             Xn[s, s] = X
    #             X = Xn

    #         X.tofile(fmt.format(*a) + '.dat')

    X= profile
    
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
    
    X.tofile('ed' + '.dat')

    dt = dxy * dt_ratio
    
    import time
    start = time.time()
    subprocess.run(['osu-hydro-pce', 'edec=0.01', 't0=0.1', 'teq=0.0', 'dt={}'.format(dt), 'dxy={}'.format(dxy), 'nls={}'.format(ls), 'tfinal=12.0', 'VisSlope = 0', 'VisHRG = 0.00', 'VisBulkMax = 0', 'VisMin = 0.0', 'InitialURead = 0'],  #'VisSlope = 0', 'VisHRG = 0.08', 'VisBulkMax = 0', 'VisMin = 0.08' 
                 check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    #print(subprocess.getoutput(['osu-hydro-pce', 'edec=0.01', 't0=0.1', 'teq=0.0', 'dt={}'.format(dt), 'dxy={}'.format(dxy), 'nls={}'.format(ls), 'tfinal=10.0', 'VisSlope = 0', 'VisHRG = 0.00', 'VisBulkMax = 0', 'VisMin = 0.0', 'InitialURead = 0']))
    end = time.time()
    print(end-start)
    surface = np.fromfile('surface.dat', dtype='f8').reshape(-1, 7)
  #  surface.tofile('surfaceAlt.dat')
    
    return dict(
    zip(['x','v'], np.hsplit(surface, [3,5])),
    #pi=dict(zip(['xx','xy','yy'],surface.T[11:14])),
    #Pi=surface.T[15],
    Temp = surface.T[5],
    #ed = surface.T[17],
    Tprop=surface.T[6])#,
   # sd = surface.T[19],
   # intersect = surface.T[20])

def run_single_event(ic):
    """
    Run the initial condition event contained in HDF5 dataset object `ic`
    and save observables to `results`.

    """    
    fs = 1#freestream.FreeStreamer(ic, grid_max, 0.1)
    # rmax = math.sqrt((run_hydro(fs, event_size=27, coarse=3)['x'][:,1:3]**2).sum(axis=1).max())
    # print(rmax)
    results = run_hydro(fs, event_size=15, profile=ic)
    return results

""" Read in profile """
profile = np.load('profiles/profilex15n13p0grid300central20.npy')

results = run_single_event(15*profile)
print("1")
#surface = frzout.Surface(**results, ymax=2)   
tvals = results['x'][:,0]
xvals = results['x'][:,1]
yvals = results['x'][:,2]
rvals = np.sqrt(xvals**2 + yvals**2)

""" Calculating fugacity at each time """    
def fugacity(t):
    return 1 - np.exp((0.1 - t)/5.0)

fugvals = fugacity(results['Tprop'])


""" tempvals gives temperature in MeV, evals gives energy density in GeV/fm^3 """
tempvals = 1000*results['Temp']
#evals = results['ed']
#intersect = results["intersect"]
Tpropvals = results['Tprop']
vvals = (results['v']**2).sum(axis=1)

plt.rcdefaults()
plt.style.use(['seaborn-darkgrid', 'seaborn-deep', 'seaborn-notebook'])
plt.rcParams.update({
    'lines.linewidth': 1.5,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Lato'],
    'mathtext.fontset': 'custom',
    'mathtext.default': 'it',
    'mathtext.rm': 'sans',
    'mathtext.it': 'sans:italic:medium',
    'mathtext.cal': 'sans',
    'pdf.fonttype': 42,
})
plt.figure(figsize=(7,5))

""" Initial condition tests """
#ic = np.fromfile('ed.dat').reshape(grid_n + 1, grid_n + 1)
#plt.plot(ic[51,:])

# rg = np.array([])
# for i in np.linspace(-8.9, 8.9, 90):
#     for j in np.linspace(-8.9, 8.9, 90):
#         rg = np.append(rg, np.sqrt(i**2 + j**2))

# glauber = np.zeros([90,90])
# with open('GlauberMC.dat') as f:
#     glaubervals = f.readlines()[1:]
#     f.close()
    
# for vals in glaubervals:
#     x = float(vals[5:15])
#     y = float(vals[15:25])
#     e = float(vals[25:])
#     if (x not in np.linspace(-8.9, 8.9, 90).round(2)):
#         print("x = ", x)
#         break
#     elif not (y in np.linspace(-8.9, 8.9, 90).round(2)):
#         print("y = ", y)
#         break
#     else:
#         xpos = int(round((x+8.9)/0.2))
#         ypos = int(round((y+8.9)/0.2))
#         glauber[xpos, ypos] = e

""" 2D scatter plot of temperature (tempvals) over space and time: """
plt.xlabel('x (fm)')
plt.ylabel(r'$\tau$ (fm/c)')
print("2")
plt.scatter(rvals[::50], tvals[::50], c = vvals[::50], cmap=plt.cm.jet)#, vmin = 0, vmax = 400)
#plt.ylim(0,16)
plt.colorbar(label=r'$v$', extend='both')




""" Plotting contours """
r = np.linspace(0, 15, 200)
t = np.linspace(0, 12, 200)
plt.title(r'$T_{eq} = 0 fm/c$')
            

# temp = griddata((rvals[::50],tvals[::50]),fugvals[::50],(r[None,:],t[:,None]),method='nearest')
# cs = plt.contour(r, t, temp, levels=[0.25, 0.5, 0.75, 0.9], colors='k', linewidths = 0.5, extend='both')
# plt.clabel(cs, inline=0, fontsize=10)


""" 2D scatter plot of fugacity (fugvals) over space and time: """
# plt.xlabel('x (fm)')
# plt.ylabel(r'$\tau$ (fm/c)')
# plt.scatter(rvals, tvals, c = fugvals, cmap=plt.cm.jet, vmin = 0, vmax = 1)
# plt.colorbar(label=r'$\lambda$')


# """ Plotting contours """
# r = np.linspace(0, 15, 300)
# t = np.linspace(0, 15, 300)
# plt.title(r'$T_{eq} = 5 fm/c$')

# fug = griddata((rvals,tvals),fugvals,(r[None,:],t[:,None]),method='linear')
# cs = plt.contour(r, t, fug, levels=[0.25, 0.5, 0.75, 0.9], colors='k', linewidths = 0.5)
# plt.clabel(cs, inline=0, fontsize=10)


# inter = griddata((rvals,tvals),intersect,(r[None,:],t[:,None]),method='linear')
# cs = plt.contour(r, t, inter, levels=[1], colors='w', linewidths = 0.5)
# plt.clabel(cs, inline=0, fontsize=10)


# plt.scatter(rvals, tvals, c = ((results['v']**2).sum(axis=1))**(1/2), cmap=plt.cm.jet, vmin = 0, vmax = 1)
# plt.colorbar(label=r'$|v|$')

# temp = griddata((rvals,tvals),tempvals,(r[None,:],t[:,None]),method='linear')
# cs2 = plt.contour(r, t, temp, levels=[155], colors='w', linewidths = 0.5)
# plt.clabel(cs2, inline=0, fontsize=10)


""" 1D plot of initial energy density over x: """
### Saved to plots/InitialEnergyDensity.png
#plt.xlabel('x (fm)')
#plt.ylabel(r'$\epsilon (GeV/fm^3)$')
#plt.scatter(rvals, evals, s = 1.0)
    
""" Plot of energy density in central cell: """
### Saved to plots/CentralCell.png
# plt.xlabel(r'$\tau$ (fm/c)'$)
# plt.ylabel(r'$\epsilon (GeV/fm^3)$')
# rvals_original = rvals
# rvals = rvals[rvals<2*grid_step]
# tvals = tvals[rvals_original<2*grid_step]
# evals = evals[rvals_original<2*grid_step]
# plt.plot(tvals, evals)


""" Determining and plotting total entropy per space-time rapidity """
# tspace = np.linspace(0.1, 12.1, 121)
# indices = []

# for t in tspace:
#     newindex = int(np.argmax(results['x'][:,0]>t))
#     if (newindex == 0 and t > tspace[0]):
#         break
#     indices.append(newindex)

# S = np.ones(len(indices))
# sd = results['sd']
# gamma = (1-(results['v']**2).sum(axis=1))**(-1/2)
# dtspace = tspace[1]-tspace[0]
# r = np.linspace(0, 15, 100)
# t = np.linspace(0, 15, 100)
# #entropy = griddata((rvals,tvals),sd*gamma*tvals,(r[None,:],t[:,None]),method='linear')

# for i in range(1,len(indices)):
#     S[i-1] = (0.001*gamma[indices[i-1]:indices[i]]*tvals[indices[i-1]:indices[i]]*sd[indices[i-1]:indices[i]]/(indices[i]-indices[i-1])).sum()
# S[len(indices)-1] = (0.001*gamma[indices[-1]:]*tvals[indices[-1]:]*sd[indices[-1]:]/(tvals.size - indices[-1])).sum()
    
# plt.xlabel(r'$\tau$ (fm/c)')
# plt.ylabel(r'$10^{-3} dS/d\eta$')
# plt.plot(tspace[:len(indices)], S)