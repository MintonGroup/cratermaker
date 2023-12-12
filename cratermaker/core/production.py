import numpy as np
from numpy.random import Generator
from scipy.optimize import fsolve
from scipy import integrate
from cratermaker.core.target import Target
from cratermaker.utils.general_utils import float_like
from numpy.typing import ArrayLike
from typing import Union, Sequence, Tuple, Callable, Any

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
    

class NeukumProduction(Production):
    
    def __init__(self, model: str = "Moon"):
        
        valid_target_name = ["Moon", "Mars", "Projectile"]
        
        if isinstance(model, str):
            if model in valid_target_name:
                self.model = model
            else:
                raise ValueError(f"Invalid model {model}. Must be one of {valid_target_name}")
        else:
            raise ValueError("model must be a string")
        
        # The Neukum production function for the Moon and Mars
        # Lunar PF from: Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to 
        #    the Lunar Reference System. Space Science Reviews 96, 55–86. https://doi.org/10.1023/A:1011989004263
        # Mars PF from: Ivanov, B.A., 2001. Mars/Moon Cratering Rate Ratio Estimates. Space Science Reviews 96, 87–104.
        #    https://doi.org/10.1023/A:1011941121102
        # Projectile PF from: Ivanov, B.A., Neukum, G., Wagner, R., 2001. Size-Frequency Distributions of Planetary Impact Craters 
        #    and Asteroids, in: Collisional Processes in the Solar System. Springer Netherlands, Dordrecht, pp. 1–34. 
        #    https://doi.org/10.1007/978-94-010-0712-2_1

        sfd_coef = {
                "Moon" : np.array(
                    [
                        -3.0876,
                        -3.557528,
                        +0.781027,
                        +1.021521,
                        -0.156012,
                        -0.444058,
                        +0.019977,
                        +0.086850,
                        -0.005874,
                        -0.006809,
                        +8.25e-4, 
                        +5.54e-5
                    ]),
                "Mars" : np.array(
                    [
                        -3.384, 
                        -3.197,
                        +1.257,
                        +0.7915,
                        -0.4861,
                        -0.3630,
                        +0.1016,
                        +6.756e-2,
                        -1.181e-2,
                        -4.753e-3,
                        +6.233e-4,
                        +5.805e-5
                    ]),
                "Projectile" : np.array(
                    [
                        0,
                        +1.375458,
                        +1.272521e-1,
                        -1.282166,
                        -3.074558e-1,
                        +4.149280e-1,
                        +1.910668e-1,
                        -4.260980e-2,
                        -3.976305e-2,
                        -3.180179e-3,
                        +2.799369e-3,
                        +6.892223e-4,
                        +2.614385e-6,
                        -1.416178e-5,
                        -1.191124e-6
                    ]
                )
            }
        sfd_range = {
                "Moon" : np.array([0.01,1000]),
                "Mars" : np.array([0.015,362]),
                "Projectile" : np.array([0.0001, 200.0]) # Estimated based on Fig. 16 of Ivanov et al. (2001)
            }
       
        #Exponential time constant 
        self.tau = 6.93
        Nexp = 5.44e-14
        Nexp_proj = 1.0 

        time_coef = {
                "Moon" : Nexp,
                "Mars" : Nexp * 10**(sfd_coef.get("Mars")[0]) / 10**(sfd_coef.get("Moon")[0]),
                "Projectile": Nexp_proj * Nexp_proj * 10**(sfd_coef.get("Projectile")[0]) / 10**(sfd_coef.get("Moon")[0]),
            }   
        
        self.valid_time = (-4500.0,None)  # Range over which the production function is valid
        
        self.sfd_coef = sfd_coef[self.model]
        self.time_coef = time_coef[self.model]
        self.sfd_range = sfd_range[self.model]
       
        
    def function(self,
             diameter: float_like | Sequence[float_like] | ArrayLike,
             time_range: Tuple[float_like, float_like] = (-1000.0, 0.0),
             check_valid_time: bool=True
             ) -> Union[float_like, ArrayLike]:
        """
        Return the cumulative size-frequency distribution of craters over a given time range and crater diameters.

        Parameters
        ----------
        diameter : float_like or numpy array
            Crater diameter(s) in units of meters to compute corresponding cumulative number density value.
        time_range : Tuple of 2 values, default=(-1000.0,0.0)
            time range relative to the present day to compute cumulative SFD in units of My. Defaults to the last 1000 My
        check_valid_time : bool, optional (default=True)
            If True, return NaN for time values outside the valid time range

        Returns
        -------
        float_like or numpy array of float_like
            The cumulative number of craters per square meter greater than the input diameter that would be expected to form on a 
            surface over the given time range.
        """            

        return self.size_frequency_distribution(diameter) * (self.chronology(time_range[1],check_valid_time) - self.chronology(time_range[0],check_valid_time))
    
     
    def chronology(self,
             time: float_like | Sequence[float_like] | ArrayLike,
             check_valid_time: bool=True
             ) -> Union[float_like, ArrayLike]:
        """
        Return the relative number of craters produced over a given time range.

        Parameters
        ----------
        time : float_like or 
            time to compute the relative number of craters in units of My.
        check_valid_time : bool, optional (default=True)
            If True, return NaN for time values outside the valid time range

        Returns
        -------
        float_like or numpy array of float_like
            The cumulative number of craters per square meter greater than the input diameter that would be expected to form on a 
            surface over the given time range.
            
        Notes
        ----- 
        The CSFD is computed using the model of Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86 for the Moon and Mars, with 
        minor changes. Notably, there is a typo in the chronology function (Eq. 5) of the original paper. The linear term in the paper
        is given as 8.38e-4. The value should be 10**(a0), and therefore the number given in the paper is based on the "Old"
        coefficients from Neukum (1983). The correct value is 10**(-3.0876) = 8.17e-4. We compute the value from the coefficients 
        in our implementation of the chronology function.
        """     
        if np.any(time) < 0.0:
            raise ValueError("time must be positive")
        
        time_Gy = -np.array(time) * 1e-3  # Convert time range from My to Gy ago for internal functions
        
        def _N1(time: float_like | Sequence[float_like] | ArrayLike,
                check_valid_time:bool=True
                ) -> Union[float_like, ArrayLike]:
            """
            Return the N(1) value as a function of time. 

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
            N1 = self.time_coef * (np.exp(self.tau * time) - 1.0) + 10 ** (self.sfd_coef[0]) * time
            if check_valid_time:
                if self.valid_time[0] is not None:
                    max_time = -self.valid_time[0] * 1e-3
                    N1 = np.where(time <= max_time, N1, np.nan)
                if self.valid_time[1] is not None:
                    min_time = -self.valid_time[1] * 1e-3
                    N1 = np.where(time >= min_time, N1, np.nan) 
            return N1.item() if np.isscalar(time) else N1
        
        N1_reference = self._CSFD(1.0) 
        N1_values = _N1(time_Gy,check_valid_time)
        N1_values /= N1_reference
        
        return  N1_values 
    
        
    def size_frequency_distribution(self,diameter: float_like | Sequence[float_like] | ArrayLike,) -> Union[float_like, ArrayLike]:

        if np.any(diameter < 0.0):
            raise ValueError("diameter must be greater than or equal to 0.0")
              
        Dkm = diameter * 1e-3 # Convert m to km for internal functions

        if self.model == "Projectile":
            Ncumulative = R_to_CSFD(R=self._CSFD, D=Dkm) 
        else:
            Ncumulative = self._CSFD(Dkm) 
            
        return Ncumulative * 1e-6 # convert from km^-2 to m^-2    
    
        
    def _CSFD(self, Dkm: float_like | Sequence[float_like] | ArrayLike) -> Union[float_like, ArrayLike]:
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
        
        def _extrapolate_sfd(side: str = "lo") -> Union[float_like, ArrayLike]:
            """
            Return the exponent, p, and and proportionality constant, A, for  the extrapolated 
            CSFD in the form N(D) = A * D**-p. 

            Parameters
            ----------
            side : str
                The side of the range to extrapolate. Valid values are "lo" and "hi"

            Returns
            -------
            A, p
                
            """    
            if side == "lo":
                idx = 0
            elif side == "hi":
                idx = 1
            else:
                raise ValueError("side must be 'lo' or 'hi'")
            p = _dNdD(self.sfd_range[idx])
            A = self._CSFD(self.sfd_range[idx])
            return A, p


        def _dNdD(Dkm: float_like | Sequence[float_like] | ArrayLike) -> Union[float_like, ArrayLike]:
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
                    _extrapolate_sfd(side="lo")
                    return A * (p / Dkm) * (Dkm / self.sfd_range[0]) ** p 
                elif Dkm > self.sfd_range[1]:
                    _extrapolate_sfd(side="hi")
                    return A * (p / Dkm) * (Dkm / self.sfd_range[0]) ** p 
                else:
                    return sum(co * np.log10(Dkm) ** i for i, co in enumerate(dcoef))
            
            return _dNdD_scalar(Dkm) if np.isscalar(Dkm) else np.vectorize(_dNdD_scalar)(Dkm)

        
        def _CSFD_scalar(Dkm):
            if Dkm < self.sfd_range[0]:
                A, p = _extrapolate_sfd(side="lo")
                return A * (Dkm / self.sfd_range[0]) ** p
            elif Dkm > self.sfd_range[1]:
                A, p = _extrapolate_sfd(side="hi")
                return A * (Dkm / self.sfd_range[1]) ** p
            else:
                logCSFD = sum(co * np.log10(Dkm) ** i for i, co in enumerate(self.sfd_coef))
                return 10 ** logCSFD
       
        return _CSFD_scalar(Dkm) if np.isscalar(Dkm) else np.vectorize(_CSFD_scalar)(Dkm)


    # def _time_to_Nrel(self,
    #                time: float_like | Sequence[float_like] | ArrayLike, 
    #                check_valid_time:bool=True
    #                )-> Union[float_like, ArrayLike]:
    #     """
    #     Return the number density of craters at a given time relative to time = 1 Gy ago.

    #     Parameters
    #     ----------
    #     time : numpy array
    #         Time in units of 
    #     check_valid_time : bool, optional (default=True)
    #         If True, return NaN for time values outside the valid time range
            
    #     Returns
    #     -------
    #     float_like or numpy array
    #        Number density of craters at the given time 
    #     """
        
    #     return self._N1(time,check_valid_time) / self._CSFD(1.0)


    # def _Nrel_to_time(self,
    #               Nrel: float_like | Sequence[float_like] | ArrayLike,
    #               check_valid_time:bool=True
    #               )-> Union[float_like, ArrayLike]:
    #     """
    #     Return the time in  for the given number density of craters relative to that at 1 Gy ago.
    #     This is the inverse of _time_to_Nrel.

    #     Parameters
    #     ----------
    #     Nrel : numpy array
    #         number density of craters relative to that at 1 Gy ago 
    #     check_valid_time : bool, optional (default=True)
    #         If True, return NaN for time values outside the valid time range

    #     Returns
    #     -------
    #     float_like or numpy array
    #         The time in Gy ago for the given relative number density of craters. 
    #     """
        
    #     def func(time,Nrel):
    #         return self._time_to_Nrel(time,check_valid_time=False) - Nrel 
        
    #     xtol = 1e-10
    #     max_guess = self.max_time * (1.0 - xtol)
    #     x0 = np.where(Nrel < max_guess, Nrel, max_guess)
    #     root_val, infodict, ier, mesg = fsolve(func=func, x0=x0, args=(Nrel), xtol=xtol, full_output=True) 
    #     if ier == 1:
    #         if check_valid_time:
    #             root_val = np.where(root_val <= self.max_time, root_val, np.nan)
    #         retval = root_val
    #     else:
    #         retval = Nrel
    #         raise ValueError(f"_Nrel_to_time failed. {mesg}")
    #     return retval.item() if np.isscalar(Nrel) else retval


def R_to_CSFD(
              R: Callable[Union[float_like, ArrayLike], Union[float_like, ArrayLike]], 
              D: Union[float_like, ArrayLike],
              Dlim: float_like = 1e6,
              *args: Any,
            ) -> Union[float_like, ArrayLike]:
    """
    Convert R values to cumulative N values for a given D using the R-plot function.

    Parameters
    ----------
    R : R = f(D) 
        A function that computes R given D.
    D : float_like or ArrayLike
        Diameters in units of km.
    Dlim : float_like
        Upper limit on the diameter over which to evaluate the integral
    *args : Any
        Additional arguments to pass to the R function

    Returns
    -------
    float or ArrayLike
        The cumulative number of craters greater than D in diameter.
    """
    def _R_to_CSFD_scalar(R, D, Dlim, *args):
        # Helper function to integrate
        def integrand(D):
            return R(D,*args) / D**3  # This is dN/dD
        
        N = 0.0
        D_i = D
        while D_i < Dlim:
            D_next = D_i * np.sqrt(2.0) 
            D_mid = (D_i + D_next) / 2  # Mid-point of the bin
            bin_width = D_next - D_i
            R_value = integrand(D_mid)
            N += R_value * bin_width
            D_i = D_next  # Move to the next bin
    
        return N
    
    return _R_to_CSFD_scalar(R, D, Dlim, *args) if np.isscalar(D) else np.vectorize(_R_to_CSFD_scalar)(R, D, Dlim, *args)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    def plot_npf_csfd():
        fig = plt.figure(1, figsize=(8, 7))
        ax = {'Moon': fig.add_subplot(121),
            'Mars': fig.add_subplot(122)}

        tvals = [0.01,1.0,4.0]
        x_min = 1e-3
        x_max = 1e4
        y_min = 1e-12
        y_max = 1e6
        nD = 1000
        Dvals = np.logspace(np.log10(x_min), np.log10(x_max), num=nD)
        for key in ax:
            production = NeukumProduction(model=key)
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
            inrange = (Dvals >= production.sfd_range[0]) & (Dvals <= production.sfd_range[1])
            lo = Dvals < production.sfd_range[0]
            hi = Dvals > production.sfd_range[1]
            for t in tvals:
                Nvals = production.function(diameter=Dvals*1e3,time_range=(-t*1e3,0.0))
                Nvals *= 1e6 # convert from m^-2 to km^-2
                ax[key].plot(Dvals[inrange], Nvals[inrange], '-', color='black', linewidth=1.0, zorder=50)
                ax[key].plot(Dvals[lo], Nvals[lo], '-.', color='orange', linewidth=2.0, zorder=50)
                ax[key].plot(Dvals[hi], Nvals[hi], '-.', color='orange', linewidth=2.0, zorder=50)
                labeli = int(0.25*nD)
                ax[key].text(Dvals[labeli],3*Nvals[labeli],f"{t:.2f} ", ha="left", va="top",rotation=-72)

        plt.tick_params(axis='y', which='minor')
        plt.tight_layout()
        plt.show()
        
    def plot_npf_N1_vs_T():
        fig = plt.figure(1, figsize=(8, 4))
        ax = fig.add_subplot(111)
        ax.set_yscale('log')
        ax.set_ylabel('$\mathregular{N(1) (km^{-2})}$')
        ax.set_xlabel('Time (Gy ago)')
        ax.set_xlim(4.5, 0)
        #ax.set_ylim(1e-1, 1e5)
        production = NeukumProduction(model="Moon")
        tvals = np.linspace(4.5, 0.0, num=1000)
        ax.plot(tvals, production.chronology(-tvals*1e3)*production.size_frequency_distribution(diameter=1000.0)*1e6, '-', color='black', linewidth=1.0, zorder=50)
        plt.tight_layout()
        plt.show() 

    
    def plot_npf_proj_rplot():
        fig = plt.figure(1, figsize=(4, 7))
        ax = fig.add_subplot(111)

        x_min = 1e-5
        x_max = 1e3
        y_min = 1e-1
        y_max = 1e2
        nD = 1000
        Dvals = np.logspace(np.log10(x_min), np.log10(x_max), num=nD)
        production = NeukumProduction(model="Projectile")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel('$\mathregular{R}$')
        ax.set_xlabel('Projectile Diameter (km)')
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
        ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
        ax.grid(True,which="minor",ls="-",lw=0.5,zorder=5)
        ax.grid(True,which="major",ls="-",lw=1,zorder=10)
        inrange = (Dvals >= production.sfd_range[0]) & (Dvals <= production.sfd_range[1])
        lo = Dvals < production.sfd_range[0]
        hi = Dvals > production.sfd_range[1]
        t = 1.0
        Nvals = production._CSFD(Dkm=Dvals)
        ax.plot(Dvals[inrange], Nvals[inrange], '-', color='black', linewidth=1.0, zorder=50)
        ax.plot(Dvals[lo], Nvals[lo], '-.', color='orange', linewidth=2.0, zorder=50)
        ax.plot(Dvals[hi], Nvals[hi], '-.', color='orange', linewidth=2.0, zorder=50)

        plt.tick_params(axis='y', which='minor')
        plt.tight_layout()
        plt.show()          

    
    def plot_npf_proj_csfd():
        fig = plt.figure(1, figsize=(4, 7))
        ax = fig.add_subplot(111)

        x_min = 1e-5
        x_max = 1e4
        y_min = 1e-7
        y_max = 1e13
        nD = 1000
        Dvals = np.logspace(np.log10(x_min), np.log10(x_max), num=nD)
        production = NeukumProduction(model="Projectile")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel('$\mathregular{N_{>D}}$')
        ax.set_xlabel('Projectile Diameter (km)')
        ax.set_xlim(x_min, x_max)
        #ax.set_ylim(y_min, y_max)
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
        ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
        ax.grid(True,which="minor",ls="-",lw=0.5,zorder=5)
        ax.grid(True,which="major",ls="-",lw=1,zorder=10)
        inrange = (Dvals >= production.sfd_range[0]) & (Dvals <= production.sfd_range[1])
        lo = Dvals < production.sfd_range[0]
        hi = Dvals > production.sfd_range[1]
        t = 1.0
        Nvals = production.function(diameter=Dvals*1e3,time_range=(-t*1e3,0.0))
        Nvals *= 1e6 # convert from m^-2 to km^-2
        ax.plot(Dvals[inrange], Nvals[inrange], '-', color='black', linewidth=1.0, zorder=50)
        ax.plot(Dvals[lo], Nvals[lo], '-.', color='orange', linewidth=2.0, zorder=50)
        ax.plot(Dvals[hi], Nvals[hi], '-.', color='orange', linewidth=2.0, zorder=50)

        plt.tick_params(axis='y', which='minor')
        plt.tight_layout()
        plt.show()    
            
    #plot_npf_csfd()
    plot_npf_N1_vs_T()
    #plot_npf_proj_rplot()
    #plot_npf_proj_csfd()
