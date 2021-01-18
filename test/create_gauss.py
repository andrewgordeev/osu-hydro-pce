import numpy as np
from scipy.optimize import curve_fit

"""
Create a Gaussian initial condition by matching a best fit to a given cylindrically symmetric TRENTO profile.
"""

profile = np.load('profiles/profilex15n13p1.npy')
grid_size = 100

def gauss(x,x0,y0,sigma):
    p = [x0, y0, sigma]
    return p[1]*np.exp(-((x-p[0])/p[2])**2)

x = np.linspace(-15, 15, grid_size)
y = x
fit, tmp = curve_fit(gauss, x, profile[:,int(grid_size/2)], p0=[1.,1.,1.])

new_grid_size = 200
x = np.linspace(-15, 15, new_grid_size)
y = x

newprofile = np.zeros((new_grid_size, new_grid_size))
for j in range(new_grid_size):
    for i in range(new_grid_size):    
        newprofile[i,j] = gauss(np.sqrt(x[i]**2 + y[j]**2), 0, fit[1], fit[2])
        
np.save('profiles/profilex15n13p1GaussAlt.npy', newprofile)