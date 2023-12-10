import numpy as np
from numpy.random import Generator
from .target import Target
from scipy import optimize

# The Neukum production function for the Moon and Mars
#Lunar PF from: Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86
#Mars PF from: Ivanov (2001) SSR v. 96 pp. 87-104
sfd_coef = {
        "NPF_Moon" : [-3.0876,-3.557528,+0.781027,+1.021521,-0.156012,-0.444058,+0.019977,+0.086850,-0.005874,-0.006809,+8.25e-4, +5.54e-5],
        "NPF_Mars" : [-3.384, -3.197,   +1.257,   +0.7915,  -0.4861,  -0.3630,  +0.1016,  +6.756e-2,-1.181e-2,-4.753e-3,+6.233e-4,+5.805e-5]
    }
sfd_range = {
    "NPF_Moon" : [0.01,1000],
    "NPF_Mars" : [0.015,362]
}
#Exponential time constant (Ga)
tau = 6.93
Nexp = 5.44e-14

time_coef = {
    "NPF_Moon" : Nexp,
    "NPF_Mars" : Nexp * 10**(sfd_coef.get("NPF_Mars")[0]) / 10**(sfd_coef.get("NPF_Moon")[0])
}

class Production():
    """
    An operations class for computing the production function for craters and impactors.


    Parameters
    ----------
    target : Target
        The target body for the impact simulation.
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used.
    """
    


    def __init__(self, 
                target: Target | str = "Moon", 
                rng: Generator | None = None):
        if not target:
            self.target = Target("Moon")
        elif isinstance(target, str):
            try:
                self.target = Target(target)
            except:
                raise ValueError(f"Invalid target name {target}")
        elif isinstance(target, Target):
            self.target = target
        else:
            raise TypeError("The 'target' argument must be a string or a Target instance")
         
        if rng is None:
            self.rng = np.random.default_rng()
        elif isinstance(rng, Generator):
            self.rng = rng
        else:
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")

        return
    


def N1(T, pfmodel):
    """
    Return the N(1) value as a function of time for a particular production function model

    Parameters
    ----------
    T : numpy array
        Time in units of Ga
    pfmodel : string
        The production function model to use
        Currently options are 'NPF_Moon' or 'NPF_Mars'

    Returns
    -------
    N1 : numpy array
        The number of craters per square kilometer greater than 1 km in diameter
    """
    T_valid_range = np.where((T <= 4.5) & (T >= 0.0), T, np.nan)
    return time_coef.get(pfmodel) * (np.exp(tau * T_valid_range) - 1.0) + 10 ** (sfd_coef.get(pfmodel)[0]) * T_valid_range

def CSFD(Dkm,planet):
    """
    Return the cumulative size-frequency distribution for a particular production function model

    Parameters
    ----------
    Dkm : numpy array
        Diameters in units of km
    pfmodel : string
        The production function model to use
        Currently options are 'NPF_Moon' or 'NPF_Mars'

    Returns
    -------
    CSFD : numpy array
        The number of craters per square kilometer greater than Dkm in diameter at T=1 Ga
    """
    logCSFD = sum(co * np.log10(Dkm) ** i for i, co in enumerate(sfd_coef.get(planet)))
    return 10**(logCSFD)

def DSFD(Dkm,planet):
    """
    Return the differential size-frequency distribution for a particular production function model

    Parameters
    ----------
    Dkm : numpy array
        Diameters in units of km
    pfmodel : string
        The production function model to use
        Currently options are 'NPF_Moon' or 'NPF_Mars'

    Returns
    -------
    DSFD : numpy array
        The differential number of craters (dN/dD) per square kilometer greater than Dkm in diameter at T = 1 Ga
    """
    dcoef = sfd_coef.get(planet)[1:]
    logDSFD = sum(co * np.log10(Dkm) ** i for i, co in enumerate(dcoef))
    return 10**(logDSFD) * CSFD(Dkm,planet) / Dkm

def Tscale(T,planet):
    """
    Return the number density of craters at time T relative to time T = 1 Ga

    Parameters
    ----------
    T : numpy array
        Time in units of Ga
    pfmodel : string
        The production function model to use
        Currently options are 'NPF_Moon' or 'NPF_Mars'

    Returns
    -------
    Tscale : numpy array
        N1(T) / CSFD(Dkm = 1.0)
    """
    return N1(T,planet) / CSFD(1.0,planet)

def pf_csfd(T,Dkm,planet):
    """
    Return the cumulative size-frequency distribution for a particular production function model as a function of Time

    Parameters
    ----------
    T : numpy array
        Time in units of Ga
    Dkm : numpy array
        Diameters in units of km
    pfmodel : string
        The production function model to use
        Currently options are 'NPF_Moon' or 'NPF_Mars'

    Returns
    -------
    pf_csfd : numpy array
        The cumulative number of craters per square kilometer greater than Dkm in diameter at time T T
    """
    D_valid_range = np.where((Dkm >= sfd_range.get(planet)[0]) & (Dkm <= sfd_range.get(planet)[1]),Dkm,np.nan)
    return CSFD(D_valid_range,planet) * Tscale(T,planet)

def pf_dsfd(T,Dkm,planet):
    """
    Return the differential size-frequency distribution for a particular production function model as a function of Time

    Parameters
    ----------
    T : numpy array
        Time in units of Ga
    Dkm : numpy array
        Diameters in units of km
    pfmodel : string
        The production function model to use
        Currently options are 'NPF_Moon' or 'NPF_Mars'

    Returns
    -------
    pf_dsfd : numpy array
        The cumulative number of craters per square kilometer greater than Dkm in diameter at time T T
    """
    D_valid_range = np.where((Dkm >= sfd_range.get(planet)[0]) & (Dkm <= sfd_range.get(planet)[1]), Dkm, np.nan)
    return DSFD(D_valid_range, planet) * Tscale(T, planet)

def T_from_scale(TS,planet):
    """
    Return the time in Ga for the given number density of craters relative to that at 1 Ga.
    This is the inverse of Tscale

    Parameters
    ----------
    TS : numpy array
        number density of craters relative to that at 1 Ga
    pfmodel : string
        The production function model to use
        Currently options are 'NPF_Moon' or 'NPF_Mars'

    Returns
    -------
    T_from_scale : numpy array
        The time in Ga
    """
    def func(S):
        return Tscale(S,planet) - TS
    return optimize.fsolve(func, np.full_like(TS,4.4),xtol=1e-10)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    print("Tests go here")
    print(f"T = 1 Ga, N(1) = {pf_csfd(1.0,1.00,'NPF_Moon')}")
    print(f"T = 4.2 Ga, N(1) = {pf_csfd(4.2,1.00,'NPF_Moon')}")
    print("Tscale test: Should return all 1s")
    Ttest = np.logspace(-4,np.log10(4.4),num=100)
    Tres = T_from_scale(Tscale(Ttest,'NPF_Mars'),'NPF_Mars')
    print(Ttest / Tres)
    #for i,t in enumerate(Ttest):
    #    print(t,Tscale(t,'Moon'),Tres[i])

    CSFDfig = plt.figure(1, figsize=(8, 7))
    ax = {'NPF_Moon': CSFDfig.add_subplot(121),
          'NPF_Mars': CSFDfig.add_subplot(122)}

    tvals = [0.01,1.0,4.0]
    x_min = 1e-3
    x_max = 1e3
    y_min = 1e-9
    y_max = 1e3
    Dvals = np.logspace(np.log10(x_min), np.log10(x_max))
    for key in ax:
        ax[key].title.set_text(key)
        ax[key].set_xscale('log')
        ax[key].set_yscale('log')
        ax[key].set_ylabel('$\mathregular{N_{>D} (km^{-2})}$')
        ax[key].set_xlabel('Diameter (km)')
        ax[key].set_xlim(x_min, x_max)
        ax[key].set_ylim(y_min, y_max)
        ax[key].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax[key].yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
        ax[key].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax[key].xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
        ax[key].grid(True,which="minor",ls="-",lw=0.5,zorder=5)
        ax[key].grid(True,which="major",ls="-",lw=1,zorder=10)
        for t in tvals:
            prod = pf_csfd(t,Dvals,key)
            ax[key].plot(Dvals, prod, '-', color='black', linewidth=1.0, zorder=50)
            labeli = 15
            ax[key].text(Dvals[labeli],prod[labeli],f"{t:.2f} Ga", ha="left", va="top",rotation=-72)

    plt.tick_params(axis='y', which='minor')
    plt.tight_layout()
    plt.show()
