import numpy as np
from numpy.random import Generator
from scipy.optimize import fsolve
from cratermaker.core.target import Target
from cratermaker.utils.general_utils import float_like
from typing import Union, Sequence, Tuple

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
       
        
        #Exponential time constant ()
        self.tau = 6.93
        Nexp = 5.44e-14

        time_coef = {
            "Moon" : Nexp,
            "Mars" : Nexp * 10**(sfd_coef.get("Mars")[0]) / 10**(sfd_coef.get("Moon")[0])
        }   
        
        self.max_time = 4.5  # Maximum time in Gy ago for which the production function is valid
        
        self.sfd_coef = sfd_coef[self.target_name]
        self.time_coef = time_coef[self.target_name]
        self.sfd_range = sfd_range[self.target_name]
       
        
    def csfd(self,
             diameter: float_like | Sequence[float_like] | np.ndarray,
             time_range: Tuple[float_like, float_like] = (-1000.0, 0.0),
             check_valid_time: bool=True
             ) -> Union[float_like, np.ndarray]:
        """
        Return the cumulative size-frequency distribution form of the production function over a given time range

        Parameters
        ----------
        diameter : float_like or numpy array
            Crater diameter(s) to compute corresponding cumulative N values in units of m
        time_range : Tuple of 2 values, default=(-1000.0,0.0)
            time range relative to the present day to compute cumulative SFD in units of My. Defaults to the last 1 
        check_valid_time : bool, optional (default=True)
            If True, return NaN for time values outside the valid time range

        Returns
        -------
        float_like or numpy array of float_like
            The cumulative number of craters per square meter greater than the input diameter that would be expected to form on a 
            surface over the given time range.
            
        Notes
        ----- 
        """
       
        if time_range[0] > time_range[1]:
            raise ValueError("time_range[0] must be less than time_range[1]")
         
        Dkm = diameter * 1e-3           # Convert m to km for internal functions
        time_Gy = -time_range * 1e-3    # Convert time range from My to Gy ago for internal functions
        
        N0 = self.CSFD(Dkm) * self.time_to_N(time_Gy[0],check_valid_time) * 1e6  
        N1 = self.CSFD(Dkm) * self.time_to_N(time_Gy[1],check_valid_time) * 1e6 
         
        return (N0 - N1)*1e6 # Convert from km^-2 to m^-2
    
    def csfd_to_time(self,
                     diameter: float_like | Sequence[float_like] | np.ndarray,
                     csfd: float_like | Sequence[float_like] | np.ndarray,
                     check_valid_time: bool=True
                     ) -> Union[float_like, np.ndarray]:
        """
        Return the time in the past for the given cumulative size-frequency distribution of craters as a function of diameter.
        
    
        """


    def _N1(self,
            time: float_like | Sequence[float_like] | np.ndarray,
            check_valid_time:bool=True
            ) -> Union[float_like, np.ndarray]:
        """
        Return the N(1) value as a function of time for a particular production function model

        Parameters
        ----------
        time : float_like or numpy array
            Time ago in units of 
        check_valid_time : bool, optional (default=True)
            If True, return NaN for time values outside the valid time range        

        Returns
        -------
        float_like or numpy array
            The number of craters per square kilometer greater than 1 km in diameter
        """
        retval =self.time_coef * (np.exp(self.tau * time) - 1.0) + 10 ** (self.sfd_coef[0]) * time
        if check_valid_time:
            retval = np.where(time <= self.max_time, retval, np.nan)
        return retval.item() if np.isscalar(time) else retval
    

    def _CSFD(self,
              Dkm: float_like | Sequence[float_like] | np.ndarray
              ) -> Union[float_like, np.ndarray]:
        """
        Return the cumulative size-frequency distribution at the reference time of 1 Gy ago. For diameter values outside 
        the range of the NPF, the CSFD is extrapolated using a power law.

        Parameters
        ----------
        Dkm : float_like or numpy array
            Diameters in units of km

        Returns
        -------
        float_like or numpy array
            The number of craters per square kilometer greater than Dkm in diameter at time=1 Gy ago.
        """
        def _CSFD_scalar(Dkm):
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
        return _CSFD_scalar(Dkm) if np.isscalar(Dkm) else np.vectorize(_CSFD_scalar)(Dkm)


    def _dNdD(self, 
              Dkm: float_like | Sequence[float_like] | np.ndarray
              ) -> Union[float_like, np.ndarray]:
        """
        Return the derivative of the cumulative size-frequency distribution as a function of diameter. For diameter values outside 
        the range of the NPF, the derivative is extrapolated using a power law.

        Parameters
        ----------
        Dkm : float_like or numpy array
            Diameters in units of km

        Returns
        -------
        float_like or numpy array
            The differential number of craters (dN/dD) per square kilometer greater than Dkm in diameter at time = 1 Gy ago.
        """        
       
        def _dNdD_scalar(Dkm): 
            dcoef = self.sfd_coef[1:]
            if Dkm < self.sfd_range[0]:
                D = self.sfd_range[0]
            elif Dkm > self.sfd_range[1]:
                D = self.sfd_range[1]
            else:
                D = Dkm
            return sum(co * np.log10(D) ** i for i, co in enumerate(dcoef))
        
        return _dNdD_scalar(Dkm) if np.isscalar(Dkm) else np.vectorize(_dNdD_scalar)(Dkm)


    def _DSFD(self, 
              Dkm: float_like | Sequence[float_like] | np.ndarray
              ) -> Union[float_like, np.ndarray]:
        """
        Return the differential size-frequency distribution of craters as a function of diameter. 

        Parameters
        ----------
        Dkm : numpy array
            Diameters in units of km

        Returns
        -------
        DSFD : numpy array
            The differential number of craters (dN/dD) per square kilometer greater than Dkm in diameter at time = 1 
        """

        return self._dNdD(Dkm) * self._CSFD(Dkm) / Dkm 


    def _time_to_Nrel(self,
                   time: float_like | Sequence[float_like] | np.ndarray, 
                   check_valid_time:bool=True
                   )-> Union[float_like, np.ndarray]:
        """
        Return the number density of craters at a given time relative to time = 1 Gy ago.

        Parameters
        ----------
        time : numpy array
            Time in units of 
        check_valid_time : bool, optional (default=True)
            If True, return NaN for time values outside the valid time range
            
        Returns
        -------
        float_like or numpy array
           Number density of craters at the given time 
        """
        
        return self._N1(time,check_valid_time) / self._CSFD(1.0)


    def _Nrel_to_time(self,
                  Nrel: float_like | Sequence[float_like] | np.ndarray,
                  check_valid_time:bool=True
                  )-> Union[float_like, np.ndarray]:
        """
        Return the time in  for the given number density of craters relative to that at 1 Gy ago.
        This is the inverse of _time_to_Nrel.

        Parameters
        ----------
        Nrel : numpy array
            number density of craters relative to that at 1 Gy ago 
        check_valid_time : bool, optional (default=True)
            If True, return NaN for time values outside the valid time range

        Returns
        -------
        float_like or numpy array
            The time in Gy ago for the given relative number density of craters. 
        """
        
        def func(time,Nrel):
            return self._time_to_Nrel(time,check_valid_time=False) - Nrel 
        
        xtol = 1e-10
        max_guess = self.max_time * (1.0 - xtol)
        x0 = np.where(Nrel < max_guess, Nrel, max_guess)
        root_val, infodict, ier, mesg = fsolve(func=func, x0=x0, args=(Nrel), xtol=xtol, full_output=True) 
        if ier == 1:
            if check_valid_time:
                root_val = np.where(root_val <= self.max_time, root_val, np.nan)
            retval = root_val
        else:
            retval = Nrel
            raise ValueError(f"_Nrel_to_time failed. {mesg}")
        return retval.item() if np.isscalar(Nrel) else retval


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    def plot_npf_csfd():
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
                ax[key].text(Dvals[labeli],prod[labeli],f"{t:.2f} ", ha="left", va="top",rotation=-72)

        plt.tick_params(axis='y', which='minor')
        plt.tight_layout()
        plt.show()
        
    def plot_npf_N1_vs_T():
        fig = plt.figure(1, figsize=(8, 4))
        ax = fig.add_subplot(111)
        ax.set_yscale('log')
        ax.set_ylabel('$\mathregular{N_{1} (km^{-2})}$')
        ax.set_xlabel('Time ()')
        ax.set_xlim(4.5, 0)
        #ax.set_ylim(1e-1, 1e5)
        npf = NeukumProductionFunction(target_name="Moon")
        tvals = np.linspace(4.5, 0.0, num=1000)
        v_N1 = np.vectorize(npf.N1)
        N1 = v_N1(tvals)
        ax.plot(tvals, N1, '-', color='black', linewidth=1.0, zorder=50)
        plt.tight_layout()
        plt.show() 
        
    plot_npf_csfd()
    plot_npf_N1_vs_T()
