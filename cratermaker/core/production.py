import numpy as np
from numpy.random import Generator
from scipy import optimize
from cratermaker.core.target import Target

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
    

class NeukumProductionFunction():
    
    def __init__(self, target_name: str = "Moon"):
        
        valid_target_name = ["Moon", "Mars"]
        
        if isinstance(target_name, str):
            if target_name in valid_target_name:
                self.target_name = target_name
            else:
                raise ValueError(f"Invalid target_name {target_name}. Must be one of {valid_target_name}")
        else:
            raise ValueError("target_name must be a string")
        
        # The Neukum production function for the Moon and Mars
        #Lunar PF from: Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86
        #Mars PF from: Ivanov (2001) SSR v. 96 pp. 87-104
        sfd_coef = {
                "Moon" : [-3.0876,-3.557528,+0.781027,+1.021521,-0.156012,-0.444058,+0.019977,+0.086850,-0.005874,-0.006809,+8.25e-4, +5.54e-5],
                "Mars" : [-3.384, -3.197,   +1.257,   +0.7915,  -0.4861,  -0.3630,  +0.1016,  +6.756e-2,-1.181e-2,-4.753e-3,+6.233e-4,+5.805e-5]
            }
        sfd_range = {
            "Moon" : [0.01,1000],
            "Mars" : [0.015,362]
        }
        
        time_range = {
            "Moon" : [0.0,4.5],
            "Mars" : [0.0,4.0]
        }
        
        #Exponential time constant (Ga)
        self.tau = 6.93
        Nexp = 5.44e-14

        time_coef = {
            "Moon" : Nexp,
            "Mars" : Nexp * 10**(sfd_coef.get("Mars")[0]) / 10**(sfd_coef.get("Moon")[0])
        }   
        
        self.sfd_coef = sfd_coef[self.target_name]
        self.time_coef = time_coef[self.target_name]
        self.sfd_range = sfd_range[self.target_name]
        self.time_range = time_range[self.target_name]
        
        
    def N1(self,T):
        """
        Return the N(1) value as a function of time for a particular production function model

        Parameters
        ----------
        T : numpy array
            Time in units of Ga

        Returns
        -------
        N1 : numpy array
            The number of craters per square kilometer greater than 1 km in diameter
        """
        if T > self.time_range[0] and T <= self.time_range[1]:
            return self.time_coef * (np.exp(self.tau * T) - 1.0) + 10 ** (self.sfd_coef[0])
        else:
            return np.nan


    def CSFD(self, Dkm):
        """
        Return the cumulative size-frequency distribution for a particular production function model
        over an extended range using power law extrapolation.

        Parameters
        ----------
        Dkm : numpy array
            Diameters in units of km

        Returns
        -------
        CSFD : numpy array
            The number of craters per square kilometer greater than Dkm in diameter at T=1 Ga
        """
        if Dkm < self.sfd_range[0]:
            slope = self.dNdD(self.sfd_range[0])
            A = self.CSFD(self.sfd_range[0])
            return A * (Dkm / self.sfd_range[0]) ** slope
        elif Dkm > self.sfd_range[1]:
            slope = self.dNdD(self.sfd_range[1])
            A = self.CSFD(self.sfd_range[1])
            return A * (Dkm / self.sfd_range[1]) ** slope
        else:
            logCSFD = sum(co * np.log10(Dkm) ** i for i, co in enumerate(self.sfd_coef))
            return 10 ** logCSFD

    def dNdD(self, Dkm):
        """
        Return the derivative of the cumulative size-frequency distribution as a function of diameter

        Parameters
        ----------
        Dkm : numpy array
            Diameters in units of km

        Returns
        -------
        dNdD : numpy array
            The differential number of craters (dN/dD) per square kilometer greater than Dkm in diameter at T = 1 Ga
        """        
        
        dcoef = self.sfd_coef[1:]
        if Dkm < self.sfd_range[0]:
            D = self.sfd_range[0]
        elif Dkm > self.sfd_range[1]:
            D = self.sfd_range[1]
        else:
            D = Dkm
        
        return sum(co * np.log10(Dkm) ** i for i, co in enumerate(dcoef))

    def DSFD(self, Dkm):
        """
        Return the differential size-frequency distribution for a particular production function model
        over an extended range using power law extrapolation.

        Parameters
        ----------
        Dkm : numpy array
            Diameters in units of km

        Returns
        -------
        DSFD : numpy array
            The differential number of craters (dN/dD) per square kilometer greater than Dkm in diameter at T = 1 Ga
        """

        return self.dNdD(Dkm) * self.CSFD(Dkm) / Dkm 


    def Tscale(self,T):
        """
        Return the number density of craters at time T relative to time T = 1 Ga

        Parameters
        ----------
        T : numpy array
            Time in units of Ga

        Returns
        -------
        Tscale : numpy array
            N1(T) / CSFD(Dkm = 1.0)
        """
        v_N1 = np.vectorize(self.N1)
        v_CSFD = np.vectorize(self.CSFD)
        
        return v_N1(T) / v_CSFD(1.0)


    def production_csfd(self,T,Dkm):
        """
        Return the cumulative size-frequency distribution for a particular production function model as a function of Time

        Parameters
        ----------
        T : numpy array
            Time in units of Ga
        Dkm : numpy array
            Diameters in units of km

        Returns
        -------
        production_csfd : numpy array
            The cumulative number of craters per square kilometer greater than Dkm in diameter at time T T
        """
        v_CSFD = np.vectorize(self.CSFD)
         
        return v_CSFD(Dkm) * self.Tscale(T)


    def production_dsfd(self,T,Dkm):
        """
        Return the differential size-frequency distribution for a particular production function model as a function of Time

        Parameters
        ----------
        T : numpy array
            Time in units of Ga
        Dkm : numpy array
            Diameters in units of km

        Returns
        -------
        production_dsfd : numpy array
            The cumulative number of craters per square kilometer greater than Dkm in diameter at time T T
        """
        v_DSFD = np.vectorize(self.DSFD)
        return v_DSFD(Dkm) * self.Tscale(T)


    def T_from_scale(self,TS):
        """
        Return the time in Ga for the given number density of craters relative to that at 1 Ga.
        This is the inverse of Tscale

        Parameters
        ----------
        TS : numpy array
            number density of craters relative to that at 1 Ga

        Returns
        -------
        T_from_scale : numpy array
            The time in Ga
        """
        def func(S):
            return self.Tscale(S) - TS
        return optimize.fsolve(func, np.full_like(TS,4.4),xtol=1e-10)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    print("Tests go here")
    npf = NeukumProductionFunction(target_name="Moon")
    print(f"T = 1 Ga, N(1) = {npf.production_csfd(1.0,1.00)}")
    print(f"T = 4.2 Ga, N(1) = {npf.production_csfd(4.2,1.00)}")
    print("Tscale test: Should return all 1s")
    Ttest = np.logspace(-4,np.log10(4.4),num=100)
    Tres = npf.T_from_scale(npf.Tscale(Ttest))
    print(Ttest / Tres)

    CSFDfig = plt.figure(1, figsize=(8, 7))
    ax = {'Moon': CSFDfig.add_subplot(121),
          'Mars': CSFDfig.add_subplot(122)}

    tvals = [0.01,1.0,4.0]
    x_min = 1e-5
    x_max = 1e4
    y_min = 1e-12
    y_max = 1e7
    nD = 1000
    Dvals = np.logspace(np.log10(x_min), np.log10(x_max), num=nD)
    for key in ax:
        npf = NeukumProductionFunction(target_name=key)
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
        inrange = (Dvals >= npf.sfd_range[0]) & (Dvals <= npf.sfd_range[1])
        lo = Dvals < npf.sfd_range[0]
        hi = Dvals > npf.sfd_range[1]
        for t in tvals:
            prod = npf.production_csfd(t,Dvals)
            ax[key].plot(Dvals[inrange], prod[inrange], '-', color='black', linewidth=1.0, zorder=50)
            ax[key].plot(Dvals[lo], prod[lo], '-.', color='orange', linewidth=1.0, zorder=50)
            ax[key].plot(Dvals[hi], prod[hi], '-.', color='orange', linewidth=1.0, zorder=50)
            labeli = int(0.5*nD)
            ax[key].text(Dvals[labeli],prod[labeli],f"{t:.2f} Ga", ha="left", va="top",rotation=-72)

    plt.tick_params(axis='y', which='minor')
    plt.tight_layout()
    plt.show()
