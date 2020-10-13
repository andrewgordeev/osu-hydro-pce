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

from scipy.special import expi
from scipy.interpolate import PchipInterpolator

import frzout


__doc__ = """
Generate an equation of state (EoS) from the hadron resonance gas EoS (at low
temperature) and the SU(3) lattice EoS (at high temperature) by connecting
their trace anomalies near the crossover range.  Print a table with columns
(e, p, s, T) or write a binary file to be read by osu-hydro.
"""


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

def _SU3(
        T, quantity='p_T4',
        Tc=0.260, anp=0.69, bnp=3.64, cnp=0.69
):
    """
            (e - 3p)/T^4 = T * d/dT(p/T^4)
    """
    t = T/Tc
    t2 = t*t

    t_alphas = 4.138 + 2*np.log(t)
    t_alphas2 = t_alphas*t_alphas
    t_alphas3 = t_alphas2*t_alphas
    alphas = 1.1424/t_alphas - 0.963 * np.log(t_alphas)/t_alphas2 + 0.4143/t_alphas3 \
    - 0.8118 * np.log(t_alphas)/t_alphas3 + 0.8118 * (np.log(t_alphas)**2)/t_alphas3

    if quantity == 'p_T4':
        # pnp found from integrating Inp below
        pnp = -anp/(2*t2) + bnp*cnp**2/2 * expi(-cnp * t) + bnp*np.exp(-cnp * t) * (cnp*t - 1)/(2*t2)
        pp = 8*np.pi**2/45 * (1 - 1.1937 * alphas + 5.3876 * alphas**(3/2) + 16.2044 * alphas**2 + 6.8392 * alphas**2 * np.log(alphas) \
                              - 45.68 * alphas**(5/2) - 36.599 * alphas**3 * np.log(alphas) + 41.896 * alphas**3 \
                              + 0.03225 * (-3526) * alphas**3)
        
        return pnp + pp
    
    elif quantity == 'e3p_T4':
        # note: with t = T/Tc,
        #   T * d/dT = T * 1/Tc d/dt = t * d/dt
        Inp = (anp + bnp*np.exp(-cnp*t))/t2
        
        t_alphas4 = t_alphas3*t_alphas
        dalphas = 2/t*(-1.1424/t_alphas2 - 0.963 * (1 - 2*np.log(t_alphas))/t_alphas3 - 0.4143 * 3/t_alphas4 \
                        - 0.8118 * (1 - 3*np.log(t_alphas))/t_alphas4 + 0.8118 * (2 - 3*np.log(t_alphas))*np.log(t_alphas)/t_alphas4)
        Ip = 8*np.pi**2/45 * t * dalphas * (-1.1937 + 5.3876 * 3/2 * alphas**(1/2) + 16.2044 * 2 * alphas \
                                            + 6.8392 * alphas * (1 + 2*np.log(alphas)) - 45.68 * 5/2 * alphas**(3/2)  \
                                            - 36.599 * alphas**2 * (1 + 3*np.log(alphas)) + 41.896 * 3 * alphas**2 \
                                            + 0.03225 * (-3526) * 3 * alphas**2)
                                                
        return Inp + Ip
    else:
        raise ValueError('unknown quantity: {}'.format(quantity))


def p_T4_lattice(T):
    """
    Lattice pressure p/T^4.

    """
    return _SU3(T, quantity='p_T4')

def p_T4_lattice_QGP(T):
    """
    Lattice pressure p/T^4.

    """
    return _hotqcd(T, quantity='p_T4')


def e3p_T4_lattice(T):
    """
    Lattice trace anomaly (e - 3p)/T^4.

    """
    return _SU3(T, quantity='e3p_T4')

def e3p_T4_lattice_QGP(T):
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


def plot(T, e3p_T4, p_T4, e_T4, args, hrg_kwargs):
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

    fig, axes_arr = plt.subplots(nrows=1, figsize=(7, 5))
    iter_axes = axes_arr#iter(axes_arr)

    @contextmanager
    def axes(title=None, ylabel=None):
        ax = iter_axes#next(iter_axes)
        yield ax
        ax.set_xlim(args.Tmin, args.Tmax)
        ax.set_xlabel('$T$ [GeV]')
        ax.axvspan(args.Ta, args.Tb, color='k', alpha=.07)
        if title is not None:
            ax.set_title(title)
        if ylabel is not None:
            ax.set_ylabel(ylabel)

    Thrg = np.linspace(T[0], args.Tb + .02, 100)
    hrg = HRGEOS(Thrg, **hrg_kwargs)
    p_T4_hrg = hrg.p_T4()
    e_T4_hrg = hrg.e_T4()
    e3p_T4_hrg = e_T4_hrg - 3*p_T4_hrg

    Tlat = np.linspace(args.Ta - .02, T[-1], 1000)
    e3p_T4_lat = e3p_T4_lattice_QGP(Tlat)
    p_T4_lat = p_T4_lattice_QGP(Tlat)
    e_T4_lat = e3p_T4_lat + 3*p_T4_lat

    with axes('Trace anomaly', '$(\epsilon - 3p)/T^4$') as ax:
       # ax.plot(Thrg, e3p_T4_hrg, **ref_line)
       # ax.plot(Tlat, e3p_T4_lat, **ref_line)
        ax.plot(T, e3p_T4)
        ax.set_ylim(0, 5)

    # with axes('Speed of sound', '$c_s^2$') as ax:
    #    # ax.plot(Thrg, hrg.cs2(), **ref_line)

    #     p = p_T4_lat*Tlat**4
    #     e = e_T4_lat*Tlat**4
    #    # ax.plot(Tlat, CubicSpline(e, p)(e, nu=1), **ref_line)

    #     p = p_T4*T**4
    #     e = e_T4*T**4
    #     ax.plot(
    #         T, CubicSpline(e, p)(e, nu=1),
    #         label='$\partial p/\partial\epsilon$'
    #     )
    #     p_spline = CubicSpline(T, p)
    #     ax.plot(
    #         T, p_spline(T, nu=1)/p_spline(T, nu=2)/T,
    #         label='$1/T\,(\partial p/\partial T)/(\partial^2p/\partial T^2)$',
    #         **comp_line
    #     )

    #     ax.set_ylim(.1, 1/3)
    #     ax.legend(loc='upper left')

    # with axes(
    #         'Other thermodynamic quantities',
    #         '$s/T^3$, $\epsilon/T^4$, $3p/T^4$'
    # ) as ax:

    #     #for T_, p_, e_ in [
    #     #        (Thrg, p_T4_hrg, e_T4_hrg),
    #     #        (Tlat, p_T4_lat, e_T4_lat),
    #     #]:
    #     #    for y in [e_ + p_, e_, 3*p_]:
    #     #        ax.plot(T_, y, **ref_line)

    #     ax.plot(T, e_T4 + p_T4, label='Entropy density $(\epsilon + p)/T$')
    #     ax.plot(T, p_spline(T, nu=1)/T**3,
    #             label='Entropy density $\partial p/\partial T$', **comp_line)
    #     ax.plot(T, e_T4, label='Energy density')
    #     ax.plot(T, 3*p_T4, label=r'Pressure $\times$ 3')

    #     ax.legend(loc='upper left')

    fig.tight_layout(pad=.2, h_pad=1.)

    if args.plot == '<show>':
        plt.show()
    else:
        fig.savefig(args.plot)


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
    
    Tc = 0.270
    tlat = Tc*np.array([0.1, 0.7, 0.74, 0.78, 0.82, 0.86, 0.9, 0.94, 0.98, 1, 1.02, 1.06, 1.10, 1.14, 1.18, 1.22, 1.26, 1.30, 1.34, 1.38, 1.42, 1.46, 1.5, 2, 2.5, 3, 3.5, 4.0, 4.5, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 80, 100, 200, 300, 400, 500, 600, 800, 1000])
    Ilat = np.array([0, .0104, .0162, .0232, .0318, .0433, .0594, .0859, .1433, 1.0008, 2.078, 2.4309, 2.4837, 2.4309, 2.3426, 2.2342, 2.1145, 1.998, 1.8867, 1.7809, 1.6810, 1.5872, 1.4995, 0.8038, 0.5057, 0.3589, 0.2736, 0.2207, 0.1855, 0.1606, 0.1266, 0.1050, 0.0903, 0.0798, 0.0720, 0.0375, 0.0265, 0.0216, 0.0191, 0.0174, 0.0154, 0.0142, 0.0112, 0.01, 0.0091, 0.0085, 0.0080, 0.0073, 0.0068])
   # plat = np.array([0.0015, 0.0023, 0.0033, 0.0046, 0.0064, 0.0087, 0.0118, 0.0164, 0.0222, 0.0571, 0.1455, 0.237, 0.325, 0.4074, 0.4837, 0.5539, 0.6181, 0.677, 0.7309, 0.7804, 0.8258, 0.8675, 1.189, 1.3319, 1.4098, 1.4582, 1.491, 1.5149, 1.533, 1.5591, 1.5768, 1.5898, 1.5998, 1.6078, 1.6444, 1.6572, 1.6641, 1.6686, 1.672, 1.6767, 1.68, 1.6887, 1.693, 1.6958, 1.6977, 1.6992, 1.7014, 1.703])
        
    e3p_T4 = PchipInterpolator(tlat,Ilat)
    delta_p_T4_spline = CubicSpline(tlat, Ilat/tlat).antiderivative()
    p_T4_0 = HRGEOS(np.array([args.Tmin]), **hrg_kwargs).p_T4()[0]

    def compute_p_T4(T):
        p_T4 = delta_p_T4_spline(T)
        p_T4 += p_T4_0
        return p_T4

    p_T4 = compute_p_T4(tlat)
    e_T4 = Ilat + 3*p_T4

    # energy density at original temperature points
    e_orig = e_T4 * tlat**4 / HBARC**3

    # compute thermodynamic quantities at evenly-spaced energy density points
    # as required by osu-hydro
    e = np.linspace(e_orig[0], e_orig[-24], args.nsteps)
    T = CubicSpline(e_orig, tlat)(e)
    Tsort = T.argsort()
    T = T[Tsort]
    # e = e[Tsort]
    p_T4 = compute_p_T4(T)
    e3p_T4 = e3p_T4(T)
    e_T4 = e3p_T4 + 3*p_T4
    
    p = compute_p_T4(T) * T**4 / HBARC**3
    s = (e + p)/T

    if args.plot:
         plot(T, e3p_T4, p_T4, e_T4, args, hrg_kwargs)
         return
     
    if args.write_bin:
        with open(args.write_bin, 'wb') as f:
            for x in [e[0], e[-1], p, s, T]:
                f.write(x.tobytes())
    else:
        # output table
        fmt = 4*'{:24.16e}'
        for row in zip(e, p, s, T):
            print(fmt.format(*row))


if __name__ == "__main__":
    main()
