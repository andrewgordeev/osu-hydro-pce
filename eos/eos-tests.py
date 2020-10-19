#!/usr/bin/env python3

import argparse
import warnings
from contextlib import contextmanager

import numpy as np
from scipy.interpolate import KroghInterpolator
try:
    from scipy.interpolate import CubicSpline
except ImportError:
    from scipy.interpolate import InterpolatedUnivariateSpline as CubicSpline
import sympy as sp
from sympy import symbols, Eq, nsolve
import frzout


def _hotqcd(
        T, quantity='p_T4',
        Tc=0.154, pid=95*np.pi**2/180, ct=3.8706, t0=0.9761,
        an=-8.7704, bn=3.9200, cn=0, dn=0.3419,
        ad=-1.2600, bd=0.8425, cd=0, dd=-0.0475
):
    """
    Compute the dimensionless pressure p/T^4 or trace anomaly (e - 3p)/T^4 for
    the HotQCD EoS (http://inspirehep.net/record/1307761).

    The pressure is parametrized in Eq. (16) and Table II.  The trace anomaly
    follows from the pressure and Eq. (5):

        (e - 3p)/T^4 = T * d/dT(p/T^4)

    """
    t = T/Tc
    t2 = t*t
    t3 = t2*t
    t4 = t2*t2

    # break parametrization into three small functions for easy differentiation
    # p/T^4 = f*g/h
    f = (1 + np.tanh(ct*(t - t0)))/2
    g = pid + an/t + bn/t2 + cn/t3 + dn/t4
    h = 1   + ad/t + bd/t2 + cd/t3 + dd/t4

    if quantity == 'p_T4':
        return f*g/h
    elif quantity == 'e3p_T4':
        t5 = t3*t2

        # derivatives of (f, g, h)
        # note: with t = T/Tc,
        #   T * d/dT = T * 1/Tc d/dt = t * d/dt
        df = ct/(2*np.cosh(ct*(t - t0))**2)
        dg = -an/t2 - 2*bn/t3 - 3*cn/t4 - 4*dn/t5
        dh = -ad/t2 - 2*bd/t3 - 3*cd/t4 - 4*dd/t5

        return t*(df*g/h + f/h*(dg - g*dh/h))
    else:
        raise ValueError('unknown quantity: {}'.format(quantity))


def p_T4_lattice(T):
    """
    Lattice pressure p/T^4.

    """
    return _hotqcd(T, quantity='p_T4')


def e3p_T4_lattice(T):
    """
    Lattice trace anomaly (e - 3p)/T^4.

    """
    return _hotqcd(T, quantity='e3p_T4')


# http://physics.nist.gov/cgi-bin/cuu/Value?hbcmevf
HBARC = 0.1973269718  # GeV fm


class HRGEOS:
    """
    Hadron resonance gas equation of state at a set of temperature points.

    """
    def __init__(self, T, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self._hrgs = [frzout.HRG(t, **kwargs) for t in T]

        self._T4 = T**4 / HBARC**3

    def _calc(self, quantity):
        return np.array([getattr(hrg, quantity)() for hrg in self._hrgs])

    def p_T4(self):
        return self._calc('pressure') / self._T4

    def e_T4(self):
        return self._calc('energy_density') / self._T4

    def e3p_T4(self):
        return (
            self._calc('energy_density') - 3*self._calc('pressure')
        ) / self._T4

    def cs2(self):
        return self._calc('cs2')


def plot(T, e3p_T4, p_T4, e_T4, T_gluon, e3p_T4_gluon, p_T4_gluon, e_T4_gluon, args, hrg_kwargs):
    import matplotlib.pyplot as plt

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

    ref_line = dict(lw=1., color='.5')
    comp_line = dict(ls='dashed', color='k', alpha=.4)

    fig, axes_arr = plt.subplots(nrows=4, figsize=(7, 20))
    iter_axes = iter(axes_arr)

    @contextmanager
    def axes(title=None, ylabel=None):
        ax = next(iter_axes)
        yield ax
        ax.set_xlim(0.1, 15)
        ax.set_xlabel(r'$\tau$' '(fm/c)')
        ax.axvspan(args.Ta, args.Tb, color='k', alpha=.07)
        if title is not None:
            ax.set_title(title)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
            
    def fugacity(tau0, tau, tau_eq):
        return (1.0 - sp.exp((tau0 - tau)/tau_eq))   
    
    def npfugacity(tau0, tau, tau_eq):
        return (1.0 - np.exp((tau0 - tau)/tau_eq))   
    
    alpha = 0.656 * 3
    beta = 0.156 * 3    
    tauf, T0 = symbols('tauf T0')
    
    def eq1(tau0, tau_eq):
        return Eq(0.156 * (tau0/tauf)**(-1/3) * (1 + alpha * fugacity(tau0, tauf, tau_eq))**(1/4) - T0, 0)
    
    def eq2(tau0, tau_eq): 
        return Eq(np.pi * 6.5**2 * tau0 * 32 * np.pi**2 / 45 * T0**3 - 17000 * (1 + alpha * fugacity(tau0, tauf, tau_eq))**(3/4) / (1 + alpha*fugacity(tau0, tauf, tau_eq) - beta * fugacity(tau0, tauf, tau_eq) * sp.log(fugacity(tau0, tauf, tau_eq))), 0)
    
    tau0 = 0.1
    tau_eq = 1
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    print(T0_val)
    print(tauf_val)
    tau_vals = np.linspace(0.1, 15, 100)
    Temp1 = T0_val * (tau0/tau_vals)**(1/3) * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(-1/4)
    
    tau_eq = 5
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    print(T0_val)
    print(tauf_val)
    Temp2 = T0_val * (tau0/tau_vals)**(1/3) * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(-1/4)
    
    tau_eq = 10
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    print(T0_val)
    print(tauf_val)
    Temp3 = T0_val * (tau0/tau_vals)**(1/3) * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(-1/4)
    
    Temp0 = T0_val * (tau0/tau_vals)**(1/3) * (1 + alpha)**(-1/4)

    with axes(r'$\tau_0 = 0.1$') as ax:
        ax.set_ylabel('T (MeV)')
        ax.plot(tau_vals, 1000*Temp0/5.06, label='Equilibrium')
        ax.plot(tau_vals, 1000*Temp1/5.06, label=r'$\tau_{eq} = 1$')
        ax.plot(tau_vals, 1000*Temp2/5.06, label=r'$\tau_{eq} = 5$')
        ax.plot(tau_vals, 1000*Temp3/5.06, label=r'$\tau_{eq} = 10$')
        ax.legend(loc = 'upper right')
        
    tau0 = 0.5
    tau_eq = 1
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    print(T0_val)
    print(tauf_val)
    tau_vals = np.linspace(0.5, 15, 100)
    Temp1 = T0_val * (tau0/tau_vals)**(1/3) * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(-1/4)
    
    tau_eq = 5
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    print(T0_val)
    print(tauf_val)
    Temp2 = T0_val * (tau0/tau_vals)**(1/3) * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(-1/4)
    
    tau_eq = 10
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    print(T0_val)
    print(tauf_val)
    Temp3 = T0_val * (tau0/tau_vals)**(1/3) * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(-1/4)
    
    Temp0 = T0_val * (tau0/tau_vals)**(1/3) * (1 + alpha)**(-1/4)

    with axes(r'$\tau_0 = 0.5$') as ax:
        ax.set_ylabel('T (MeV)')
        ax.plot(tau_vals, 1000*Temp0/5.06, label='Equilibrium')
        ax.plot(tau_vals, 1000*Temp1/5.06, label=r'$\tau_{eq} = 1$')
        ax.plot(tau_vals, 1000*Temp2/5.06, label=r'$\tau_{eq} = 5$')
        ax.plot(tau_vals, 1000*Temp3/5.06, label=r'$\tau_{eq} = 10$')
        ax.legend(loc = 'upper right')

    tau0 = 0.1
    tau_eq = 1
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    tau_vals = np.linspace(0.1, 15, 100)
    SEta1 = np.pi * 6.5**2 * tau0 * 32/45 * np.pi**2 * T0_val**3 * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq) - beta * npfugacity(tau0, tau_vals, tau_eq) * np.log(npfugacity(tau0, tau_vals, tau_eq))) / (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(3/4)
  
    tau_eq = 5
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    SEta2 = np.pi * 6.5**2 * tau0 * 32/45 * np.pi**2 * T0_val**3 * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq) - beta * npfugacity(tau0, tau_vals, tau_eq) * np.log(npfugacity(tau0, tau_vals, tau_eq))) / (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(3/4)
    
    tau_eq = 10
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    SEta3 = np.pi * 6.5**2 * tau0 * 32/45 * np.pi**2 * T0_val**3 * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq) - beta * npfugacity(tau0, tau_vals, tau_eq) * np.log(npfugacity(tau0, tau_vals, tau_eq))) / (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(3/4)

    SEta0 = np.pi * 6.5**2 * tau0 * 32/45 * np.pi**2 * T0_val**3 * np.ones(SEta3.shape) * (1 + alpha)**(1/4)

    with axes(r'$\tau_0 = 0.1$') as ax:
        ax.set_ylabel(r'$10^{-3} dS/d\eta$')
        ax.plot(tau_vals, 0.001*SEta0, label='Equilibrium')
        ax.plot(tau_vals, 0.001*SEta1, label=r'$\tau_{eq} = 1$')
        ax.plot(tau_vals, 0.001*SEta2, label=r'$\tau_{eq} = 5$')
        ax.plot(tau_vals, 0.001*SEta3, label=r'$\tau_{eq} = 10$')
        ax.legend(loc = 'lower right')

    tau0 = 0.5
    tau_eq = 1
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    tau_vals = np.linspace(0.5, 15, 100)
    SEta1 = np.pi * 6.5**2 * tau0 * 32/45 * np.pi**2 * T0_val**3 * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq) - beta * npfugacity(tau0, tau_vals, tau_eq) * np.log(npfugacity(tau0, tau_vals, tau_eq))) / (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(3/4)
  
    tau_eq = 5
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    SEta2 = np.pi * 6.5**2 * tau0 * 32/45 * np.pi**2 * T0_val**3 * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq) - beta * npfugacity(tau0, tau_vals, tau_eq) * np.log(npfugacity(tau0, tau_vals, tau_eq))) / (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(3/4)
    
    tau_eq = 10
    sol = nsolve((eq1(tau0, tau_eq),eq2(tau0, tau_eq)),(tauf,T0),(1,1))
    T0_val = sol[1]
    tauf_val = sol[0]
    SEta3 = np.pi * 6.5**2 * tau0 * 32/45 * np.pi**2 * T0_val**3 * (1 + alpha * npfugacity(tau0, tau_vals, tau_eq) - beta * npfugacity(tau0, tau_vals, tau_eq) * np.log(npfugacity(tau0, tau_vals, tau_eq))) / (1 + alpha * npfugacity(tau0, tau_vals, tau_eq))**(3/4)

    SEta0 = np.pi * 6.5**2 * tau0 * 32/45 * np.pi**2 * T0_val**3 * np.ones(SEta3.shape) * (1 + alpha)**(1/4)

    with axes(r'$\tau_0 = 0.5$') as ax:
        ax.set_ylabel(r'$10^{-3} dS/d\eta$')
        ax.plot(tau_vals, 0.001*SEta0, label='Equilibrium')
        ax.plot(tau_vals, 0.001*SEta1, label=r'$\tau_{eq} = 1$')
        ax.plot(tau_vals, 0.001*SEta2, label=r'$\tau_{eq} = 5$')
        ax.plot(tau_vals, 0.001*SEta3, label=r'$\tau_{eq} = 10$')
        ax.legend(loc = 'lower right')

    fig.tight_layout(pad=.2, h_pad=1.)


def T_points(Tmin, Tmax, n, extra_low=0, extra_high=0):
    """
    Create evenly-spaced temperature points, with optional "extra" points
    outside the range.  Total number of points is (n + extra_low + extra_high).

    """
    T = np.arange(-extra_low, n + extra_high, dtype=float)
    T *= (Tmax - Tmin)/(n - 1)
    T += Tmin
    return T

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--Tmin', type=float, default=.050,
        help='minimum temperature'
    )
    parser.add_argument(
        '--Ta', type=float, default=.165,
        help='connection range minimum temperature'
    )
    parser.add_argument(
        '--Tb', type=float, default=.200,
        help='connection range maximum temperature'
    )
    parser.add_argument(
        '--Tmax', type=float, default=.500,
        help='maximum temperature'
    )
    parser.add_argument(
        '--nsteps', type=int, default=10**5,
        help='number of energy density steps'
    )
    parser.add_argument(
        '--species', choices=['all', 'urqmd', 'id'], default='urqmd',
        help='HRG species'
    )
    parser.add_argument(
        '--res-width-off', action='store_false', dest='res_width',
        help='do not account for finite width of resonances'
    )
    parser.add_argument(
        '--plot', nargs='?', metavar='FILE', const='<show>', default=False,
        help='plot EoS instead of printing table'
    )
    parser.add_argument(
        '--write-bin', metavar='FILE',
        help='write binary file instead of printing table'
    )
    args = parser.parse_args()

    hrg_kwargs = dict(species=args.species, res_width=args.res_width)

    # split full temperature range into three parts:
    #   low-T (l):  Tmin < T < Ta  (HRG)
    #   mid-T (m):  Ta < T < Tb  (connection)
    #   high-T (h):  Tb < T < Tmax  (lattice)

    # number of extra temperature points below Tmin and above Tmax
    # (helps cubic interpolation)
    nextrapts = 2

    # compute low-T (HRG) trace anomaly
    Tl = T_points(args.Tmin, args.Ta, 200, extra_low=nextrapts)
    e3p_T4_l = HRGEOS(Tl, **hrg_kwargs).e3p_T4()

    # compute mid-T (connection) using an interpolating polynomial that
    # matches the function values and first several derivatives at the
    # connection endpoints (Ta, Tb)
    nd = 5

    # use Krogh interpolation near the endpoints to estimate derivatives
    def derivatives(f, T0, dT=.001):
        # evaluate function at Chebyshev zeros as suggested by docs
        T = T0 + dT*np.cos(np.linspace(0, np.pi, nd))
        return KroghInterpolator(T, f(T)).derivatives(T0)

    # use another Krogh interpolation for the connection polynomial
    # skip (Ta, Tb) endpoints in Tm since they are already in (Tl, Th)
    Tm = T_points(args.Ta, args.Tb, 100, extra_low=-1, extra_high=-1)
    e3p_T4_m = KroghInterpolator(
        nd*[args.Ta] + nd*[args.Tb],
        np.concatenate([
            derivatives(lambda T: HRGEOS(T, **hrg_kwargs).e3p_T4(), args.Ta),
            derivatives(e3p_T4_lattice, args.Tb)
        ])
    )(Tm)

    # compute high-T part (lattice)
    Th = T_points(args.Tb, args.Tmax, 1000, extra_high=nextrapts)
    e3p_T4_h = e3p_T4_lattice(Th)

    # join temperature ranges together
    T = np.concatenate([Tl, Tm, Th])
    e3p_T4 = np.concatenate([e3p_T4_l, e3p_T4_m, e3p_T4_h])

    # pressure is integral of trace anomaly over temperature starting from some
    # reference temperature, Eq. (12) in HotQCD paper:
    #   p/T^4(T) = p/T^4(T_0) + \int_{T_0}^T dT (e - 3p)/T^5
    delta_p_T4_spline = CubicSpline(T, e3p_T4/T).antiderivative()
    p_T4_0 = HRGEOS(T[:1], **hrg_kwargs).p_T4()[0]

    def compute_p_T4(T):
        p_T4 = delta_p_T4_spline(T)
        p_T4 += p_T4_0
        return p_T4

    p_T4 = compute_p_T4(T)
    e_T4 = e3p_T4 + 3*p_T4

    # energy density at original temperature points
    e_orig = e_T4 * T**4 / HBARC**3

    # compute thermodynamic quantities at evenly-spaced energy density points
    # as required by osu-hydro
    e = np.linspace(e_orig[nextrapts], e_orig[-nextrapts - 1], args.nsteps)
    T = CubicSpline(e_orig, T)(e)
    p = compute_p_T4(T) * T**4 / HBARC**3
    s = (e + p)/T
    
    Tc = 0.270
    tlat = Tc*np.array([0.7, 0.74, 0.78, 0.82, 0.86, 0.9, 0.94, 0.98, 1, 1.02, 1.06, 1.10, 1.14, 1.18, 1.22, 1.26, 1.30, 1.34, 1.38, 1.42, 1.46, 1.5, 2, 2.5, 3, 3.5, 4.0, 4.5, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 80, 100, 200, 300, 400, 500, 600, 800, 1000])
    Ilat = np.array([.0104, .0162, .0232, .0318, .0433, .0594, .0859, .1433, 1.0008, 2.078, 2.4309, 2.4837, 2.4309, 2.3426, 2.2342, 2.1145, 1.998, 1.8867, 1.7809, 1.6810, 1.5872, 1.4995, 0.8038, 0.5057, 0.3589, 0.2736, 0.2207, 0.1855, 0.1606, 0.1266, 0.1050, 0.0903, 0.0798, 0.0720, 0.0375, 0.0265, 0.0216, 0.0191, 0.0174, 0.0154, 0.0142, 0.0112, 0.01, 0.0091, 0.0085, 0.0080, 0.0073, 0.0068])
    plat = np.array([0.0015, 0.0023, 0.0033, 0.0046, 0.0064, 0.0087, 0.0118, 0.0164, 0.0222, 0.0571, 0.1455, 0.237, 0.325, 0.4074, 0.4837, 0.5539, 0.6181, 0.677, 0.7309, 0.7804, 0.8258, 0.8675, 1.189, 1.3319, 1.4098, 1.4582, 1.491, 1.5149, 1.533, 1.5591, 1.5768, 1.5898, 1.5998, 1.6078, 1.6444, 1.6572, 1.6641, 1.6686, 1.672, 1.6767, 1.68, 1.6887, 1.693, 1.6958, 1.6977, 1.6992, 1.7014, 1.703])

    e3p_T4_gluon = CubicSpline(tlat,Ilat)
    delta_p_T4_spline_gluon = CubicSpline(tlat, Ilat/tlat).antiderivative()
    p_T4_0_gluon = HRGEOS(tlat[:1], **hrg_kwargs).p_T4()[0]

    def compute_p_T4_gluon(T):
        p_T4 = delta_p_T4_spline_gluon(T)
        p_T4 += p_T4_0_gluon
        return p_T4

    p_T4_gluon = compute_p_T4_gluon(tlat)
    e_T4_gluon = Ilat + 3*plat

    # energy density at original temperature points
    e_orig_gluon = e_T4_gluon * tlat**4 / HBARC**3

    # compute thermodynamic quantities at evenly-spaced energy density points
    # as required by osu-hydro
    e_gluon = np.linspace(e_orig_gluon[nextrapts], e_orig_gluon[-nextrapts - 1], args.nsteps)
    T_gluon = CubicSpline(e_orig_gluon, tlat)(e)
    p_gluon = compute_p_T4_gluon(T) * T_gluon**4 / HBARC**3
    s_gluon = (e_gluon + p_gluon)/T_gluon
    
    plot(T, e3p_T4, p_T4, e_T4, T_gluon, e3p_T4_gluon, p_T4_gluon, e_T4_gluon, args, hrg_kwargs)
    return

if __name__ == "__main__":
    main()
