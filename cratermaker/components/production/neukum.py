import numpy as np
from numpy.random import Generator
from cratermaker.components.production import Production
from cratermaker.utils.custom_types import FloatLike
from cratermaker.utils.general_utils import R_to_CSFD, parameter, format_large_units
from numpy.typing import ArrayLike
from collections.abc import Sequence
from typing import Any, Union

@Production.register("neukum")
class NeukumProduction(Production):
    """
    An operations class for computing the the Neukum production function for the Moon and Mars.

    Parameters
    ----------
    version : {"Moon", "Mars", "Projectile"}, optional
        The specific model to use for the production function. "Moon" and "Mars" are both crater production functions, and
        "Projectile" is a projectile function. Defaults to "Moon".
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments. 
        
    Notes
    ----- 
    The CSFD is computed using the model of Ivanov, Neukum, and Hartmann (2001) for the Moon and Mars, with 
    minor changes. Notably, there is a typo in the chronology function (Eq. 5) of the original paper. The linear term in the paper
    is given as 8.38e-4. The value should be 10^(a0), and therefore the number given in the paper is based on the "Old"
    coefficients from Neukum (1983). The correct value is 10^(-3.0876) = 8.17e-4. We compute the value from the coefficients 
    in our implementation of the chronology function.       
    
    The Mars production function is based on Ivanov (2001). The projectile production function is from Ivanov, Neukum, and Wagner (2001).
   
    .. rubric:: References

    - Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to 
      the Lunar Reference System. *Space Science Reviews*, 96, 55-86. https://doi.org/10.1023/A:1011989004263
    - Ivanov, B.A., 2001. Mars/Moon Cratering Rate Ratio Estimates. *Space Science Reviews*, 96, 87-104. https://doi.org/10.1023/A:1011941121102
    - Ivanov, B.A., Neukum, G., Wagner, R., 2001. Size-Frequency Distributions of Planetary Impact Craters 
      and Asteroids. In *Collisional Processes in the Solar System*, Springer Netherlands, Dordrecht, pp. 1-34. 
      https://doi.org/10.1007/978-94-010-0712-2_1
    - Neukum, G., 1983. Meteorite bombardment and planetary surfaces dating. Habilitation Dissertation for Faculty Membership, University of Munich.

    """
    def __init__(self, 
                 version: str | None = None,
                 rng: Generator | None = None, 
                 rng_seed: int | None = None,
                 rng_state: dict | None = None,
                 **kwargs: Any):
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs) 

        self.version = version or "Moon"

    def __repr__(self) -> str:
        base = super().__repr__()
        timelo = format_large_units(self.valid_time[0], quantity="time")
        timehi = format_large_units(self.valid_time[1], quantity="time")
        dlo = format_large_units(self.sfd_range[0], quantity="length")
        dhi = format_large_units(self.sfd_range[1], quantity="length")
        return (
            f"{base}\n"
            f"Version: {self.version}\n"
            f"Valid Time Range: {timelo} - {timehi}\n"
            f"Valid Diameter Range: {dlo} - {dhi}"
        )

    def function(self,
             diameter: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
             age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
             age_end: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
             check_valid_time: bool=True,
             **kwargs: Any,
             ) -> Union[FloatLike, ArrayLike]:
        """
        Return the cumulative size-frequency distribution of craters over a given age range and crater diameter.

        Parameters
        ----------
        diameter : FloatLike or numpy array
            Crater diameter(s) in units of meters to compute corresponding cumulative number density value.
        age : FloatLike or ArrayLike, default=1.0
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        age_end, FloatLike or ArrayLike, optional
            The ending age in units of My relative to the present, which is used to compute the cumulative SFD. The default is 0 (present day).
        check_valid_time : bool, optional (default=True)
            If True, return NaN for age values outside the valid age range

        Returns
        -------
        FloatLike or numpy array of FloatLike
            The cumulative number of craters per square meter greater than the input diameter that would be expected to form on a 
            surface over the given age range.
        """
        age, age_end = self._validate_age(age, age_end) 
        diameter, _ = self._validate_csfd(diameter=diameter)
        
        diameter_array = np.asarray(self.csfd(diameter))
        age_difference = np.asarray(self.chronology(age, check_valid_time)) - np.asarray(self.chronology(age_end, check_valid_time))

        if diameter_array.ndim > 0 and age_difference.ndim > 0:
            return diameter_array[:, None] * age_difference
        else:
            return diameter_array * age_difference
    
    def chronology(self,
             age: FloatLike | Sequence[FloatLike] | ArrayLike = 1000.0,
             check_valid_time: bool=True,
             **kwargs: Any,
             ) -> Union[FloatLike, ArrayLike]:
        """
        Returns the relative number of craters produced over a given age range. This implements the chronology function given in
        Eq. 5 of Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86, but takes in the age argument in the Cratermaker unit 
        system of My instead of Gy. 

        Parameters
        ----------
        age : FloatLike or ArrayLike, default=1000.0
            Age in the past relative to the present day to compute cumulative SFD in units of My. 
        check_valid_time : bool, optional (default=True)
            If True, return NaN for age values outside the valid age range

        Returns
        -------
        FloatLike or numpy array of FloatLike
            The number of craters relative to the amount produced in the last 1 My.
            
        """     
        time_Gy = np.array(age) * 1e-3  # Convert age range from My to Gy ago for internal functions
        
        def _N1km(age: FloatLike | Sequence[FloatLike] | ArrayLike,
                check_valid_time:bool=True
                ) -> Union[FloatLike, ArrayLike]:
            """
            Return the cumulative number of 1 km craters as a function of age in Gy. This is a direct implementation of Eq. 5 in
            Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86 (with corrected coefficient for the linear term).  

            Parameters
            ----------
            age : FloatLike or numpy array
                Time ago in units of Gy
            check_valid_time : bool, optional (default=True)
                If True, return NaN for age values outside the valid age range        

            Returns
            -------
            FloatLike or numpy array
                The number of craters per square kilometer greater than 1 km in diameter
            """
            N1 = self.Cexp * (np.exp(age/self.tau) - 1.0) + self.Clin * age
            if check_valid_time:
                if self.valid_time[0] is not None:
                    min_time = self.valid_time[0] * 1e-3
                    N1 = np.where(age >= min_time, N1, np.nan)
                if self.valid_time[1] is not None:
                    max_time = self.valid_time[1] * 1e-3
                    N1 = np.where(age <= max_time, N1, np.nan) 
            return N1.item() if np.isscalar(age) else N1
        
        N1_reference = _N1km(age=1.0e-3) 
        N1_values = _N1km(age=time_Gy,check_valid_time=check_valid_time)
        N1_values /= N1_reference
        
        return  N1_values 

    @property
    def valid_versions(self) -> list[str]:
        """
        The valid versions of the production function. 
        
        Returns
        -------
        list
            The valid versions of the production function.
        """
        return ["Moon", "Mars", "Projectile"] 

    @property
    def generator_type(self) -> str:
        """
        Get the generator type of the production function.
        
        Returns
        -------
        str
            The generator type of the production function.
        """
        return self._generator_type
    
    @generator_type.setter
    def generator_type(self, value: str) -> None:
        """
        The generator type for this production function is determined by the version. 
        """
        raise NotImplementedError("The generator type cannot be changed for this production function.")

    @parameter
    def version(self) -> str:
        """
        Get the version of the production function.
        
        Returns
        -------
        str
            The version of the production function.
        """
        return self._version
    
    @version.setter
    def version(self, value: str) -> None:
        """
        Set the version of the production function.
        
        Parameters
        ----------
        value : str
            The version of the production function.
        """
        if value.title() not in self.valid_versions:
            raise ValueError(f"Invalid version '{value}'. Valid options are {self.valid_versions}")
        self._version = value.title()

        if self._version == "Projectile":
            self._generator_type = "projectile"
        else:
            self._generator_type = "crater"

    @property
    def sfd_tables(self) -> dict[str, list[float]]:
        """
        Get the table of coefficients for the production function.
        
        Returns
        -------
        dict
            The table of coefficients for the production function. 
        """
        return {"Moon" : [
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
                    ],
                "Mars" : [
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
                    ],
                "Projectile" : [
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
            }
        
    @property
    def sfd_coef(self) -> list[float]:
        """
        Get the table of coefficients for the production function.
        
        Returns
        -------
        list
            The table of SFD coefficients for the production function. 
        """
        return self.sfd_tables[self.version]
    
    @property
    def Clin(self) -> float:
        """
        Get the linear coefficient for the production function.
        
        Returns
        -------
        float
            The linear coefficient for the production function.
        """
        return 10**(self.sfd_coef[0])
    
    @property
    def Cexp(self) -> float:
        """
        Get the exponential coefficient for the production function.
        
        Returns
        -------
        float
            The exponential coefficient for the production function.
        """
        Cexp_moon = 5.44e-14
        if self.version == "Moon":
            return Cexp_moon
        else:
            Clin_moon = 10**(self.sfd_tables["Moon"][0])
            return Cexp_moon * self.Clin / Clin_moon 
    
    @property
    def tau(self) -> float:
        """
        Get the time constant for the production function.
        
        Returns
        -------
        float
            The time constant for the production function.
        """
        return 1.0 / 6.93

    @property
    def sfd_range(self) -> tuple[float, float]:
        """
        The range of diameters over which the SFD is valid. The range is given in m. 
        
        Returns
        -------
        tuple
            The lower and upper bounds of the SFD range in m.
        """
        sfd_range = {
                "Moon" : (10.0,1000.0e3),
                "Mars" : (15.0,362.0e3),
                "Projectile" : (0.1, 200.0e3) # Estimated based on Fig. 16 of Ivanov et al. (2001)
            }

        return sfd_range[self.version]
    
    @property
    def valid_time(self) -> tuple[float, float]:
        """
        The range of ages over which the production function is valid. The range is given in My. 
        
        Returns
        -------
        tuple
            The lower and upper bounds of the valid time range in My.
        """
        return (0,4500)

    
    def csfd(self,diameter: FloatLike | ArrayLike,) -> FloatLike | ArrayLike:
        """
        Return the cumulative size frequency distribution of craters at a given age relative to age = 1 My ago per m^2.

        Parameters
        ----------
        diameter : FloatLike or ArrayLike
            units of meter 
            
        Returns
        -------
        FloatLike or numpy array
           Cumulative number density of per square meter greater than the input diameter.
        """      
        # Convert m to km for internal functions
        Dkm_lo = self.sfd_range[0] * 1e-3
        Dkm_hi = self.sfd_range[1] * 1e-3

        def _extrapolate_sfd(side: str = "lo") -> Union[FloatLike, ArrayLike]:
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
                Dkm = Dkm_lo
            elif side == "hi":
                Dkm = Dkm_hi
            else:
                raise ValueError("side must be 'lo' or 'hi'")
            p = _dNdD(Dkm)
            A = _CSFD(Dkm)
            return A, p


        def _dNdD(Dkm: FloatLike | Sequence[FloatLike] | ArrayLike) -> Union[FloatLike, ArrayLike]:
            """
            Return the derivative of the cumulative size-frequency distribution as a function of diameter. For diameter values outside 
            the range of the NPF, the derivative is extrapolated using a power law.

            Parameters
            ----------
            Dkm : FloatLike or numpy array
                diameter in units of km

            Returns
            -------
            FloatLike or numpy array
                The differential number of craters (dN/dD) per square kilometer greater than Dkm in diameter at age = 1 Gy ago.
            """        
            def _dNdD_scalar(Dkm): 
                dcoef = self.sfd_coef[1:]
                if Dkm < Dkm_lo:
                    A, p = _extrapolate_sfd(side="lo")
                    return A * (p / Dkm) * (Dkm / Dkm_lo) ** p 
                elif Dkm > Dkm_hi:
                    A, p = _extrapolate_sfd(side="hi")
                    return A * (p / Dkm) * (Dkm / Dkm_hi) ** p 
                else:
                    return sum(co * np.log10(Dkm) ** i for i, co in enumerate(dcoef))
            
            return _dNdD_scalar(Dkm) if np.isscalar(Dkm) else np.vectorize(_dNdD_scalar)(Dkm)

    
        def _CSFD(Dkm: FloatLike | Sequence[FloatLike] | ArrayLike) -> Union[FloatLike, ArrayLike]:
            """
            Return the cumulative size-frequency distribution at the reference age of 1 Gy ago. For diameter values outside 
            the range of the NPF, the CSFD is extrapolated using a power law.

            Parameters
            ----------
            Dkm : FloatLike or numpy array
                diameter in units of km

            Returns
            -------
            FloatLike or numpy array
                The number of craters per square kilometer greater than Dkm in diameter at age=1 Gy ago.
            """
            def _CSFD_scalar(Dkm):
                if Dkm < Dkm_lo:
                    A, p = _extrapolate_sfd(side="lo")
                    return A * (Dkm / Dkm_lo) ** p
                elif Dkm > Dkm_hi:
                    A, p = _extrapolate_sfd(side="hi")
                    p -= 2.0 # Steepen the upper branch of the SFD to prevent anomolously large craters from forming
                    return A * (Dkm / Dkm_hi) ** p
                else:
                    logCSFD = sum(co * np.log10(Dkm) ** i for i, co in enumerate(self.sfd_coef))
                    return 10 ** logCSFD
        
            return _CSFD_scalar(Dkm) if np.isscalar(Dkm) else np.vectorize(_CSFD_scalar)(Dkm)


        if np.any(diameter < 0.0):
            raise ValueError("diameter must be greater than or equal to 0.0")
              
        Dkm = diameter * 1e-3 # Convert m to km for internal functions

        if self.version == "Projectile":
            Ncumulative = R_to_CSFD(R=_CSFD, D=Dkm) * 2.94e-5 # This is a multiplication factor that gets the projectile CSFD to approximately match the lunar crater CSFD 
        else:
            Ncumulative = _CSFD(Dkm) 
            
        return Ncumulative * 1e-9 # convert from Gy^-1 km^-2 to My^-1 m^-2    


if __name__ == "__main__":
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from cratermaker.components.target import Target
    
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
            production = NeukumProduction(version=key)
            ax[key].title.set_text(key)
            ax[key].set_xscale('log')
            ax[key].set_yscale('log')
            ax[key].set_ylabel('$\\mathregular{N_{>D} (km^{-2})}$')
            ax[key].set_xlabel('Diameter (km)')
            ax[key].set_xlim(x_min, x_max)
            ax[key].set_ylim(y_min, y_max)
            ax[key].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
            ax[key].yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
            ax[key].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
            ax[key].xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
            ax[key].grid(True,which="minor",ls="-",lw=0.5,zorder=5)
            ax[key].grid(True,which="major",ls="-",lw=1,zorder=10)
            Dkm_lo = production.sfd_range[0] * 1e-3
            Dkm_hi = production.sfd_range[1] * 1e-3
            inrange = (Dvals >= Dkm_lo) & (Dvals <= Dkm_hi)
            lo = Dvals < Dkm_lo
            hi = Dvals > Dkm_hi
            for t in tvals:
                Nvals = production.function(diameter=Dvals*1e3,age=t*1e3)
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
        ax.set_ylabel('$\\mathregular{N(1) (km^{-2})}$')
        ax.set_xlabel('Time (Gy ago)')
        ax.set_xlim(4.5, 0)
        moon = NeukumProduction(version="Moon")
        mars = NeukumProduction(version="Mars")
        tvals = np.linspace(4.5, 0.0, num=1000)
        N1_moon = moon.function(diameter=1000.0, age=tvals*1e3)*1e6
        N1_mars = mars.function(diameter=1000.0, age=tvals*1e3)*1e6
        ax.plot(tvals, N1_moon, '-', color='dimgrey', linewidth=2.0, zorder=50, label="Moon")
        ax.plot(tvals, N1_mars, '-', color='orange', linewidth=2.0, zorder=50, label="Mars")
        ax.legend()
        plt.tight_layout()
        plt.show() 


    def plot_npf_fit():
        # Define the power-law function
        def power_law(D, N1_coef, slope):
            return N1_coef * D ** slope
       
        # Create the lunar crater and projectile production functions 
        crater_production = NeukumProduction(version="Moon")
        projectile_production = NeukumProduction(version="Projectile",impact_velocity_model = "Moon_MBA")

        # Make a population of sizes that spans the size of interest
        Dc = np.logspace(-1,2)
        Di = np.logspace(-2,1)

        # Compute the reference N values (N>D at 1 My)
        Nc = crater_production.function(diameter=Dc,age=1.0)
        Ni = projectile_production.function(diameter=Di,age=1.0)
       
        D1proj = 22.3 # Approximate diameter of crater made by a 1 m projectile using the default Cratermaker scaling relationships
         
        N1_c = crater_production.function(diameter = D1proj, age = 1.0)
        N1_p = projectile_production.function(diameter = 1.0, age = 1.0)        

        # Fit the power-law function to the crater data
        params_crater, _ = curve_fit(power_law, Dc, Nc)

        # Fit the power-law function to the projectile data
        params_projectile, _ = curve_fit(power_law, Di, Ni * N1_c / N1_p)

        # Extracting the parameters
        N1_coef_crater, slope_crater = params_crater
        N1_coef_projectile, slope_projectile = params_projectile

        # Output the results
        print("Crater Production Fit Parameters N1_coef =", N1_coef_crater, ", slope =", slope_crater)
        print("Projectile Production Fit Parameters N1_coef =", N1_coef_projectile, ", slope =", slope_projectile)
        
        # Create the fitted curve data
        fitted_curve_crater = power_law(Dc, *params_crater)
        fitted_curve_projectile = power_law(Di, *params_projectile)

        # Plotting the crater data
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.loglog(Dc, Nc, 'o', label='Original Data (Crater)')
        plt.loglog(Dc, fitted_curve_crater, '-', label='Fitted Curve (Crater)')
        plt.xlabel('Crater Diameter (m)')
        plt.ylabel('N>D')
        plt.title('Crater Production')
        plt.legend()

        # Plotting the projectile data
        plt.subplot(1, 2, 2)
        plt.loglog(Di, Ni * N1_c/N1_p, 'o', label='Original Data (Projectile)')
        plt.loglog(Di, fitted_curve_projectile, '-', label='Fitted Curve (Projectile)')
        plt.xlabel('Projectile Diameter (m)')
        plt.ylabel('N>D')
        plt.title('Projectile Production')
        plt.legend()

        # Show the plot
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
        production = NeukumProduction(version="Projectile",impact_velocity_model = "Moon_MBA")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel('$\\mathregular{N_{>D}}$')
        ax.set_xlabel('Projectile Diameter (km)')
        ax.set_xlim(x_min, x_max)
        #ax.set_ylim(y_min, y_max)
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
        ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
        ax.grid(True,which="minor",ls="-",lw=0.5,zorder=5)
        ax.grid(True,which="major",ls="-",lw=1,zorder=10)
        Dkm_lo = production.sfd_range[0] * 1e-3
        Dkm_hi = production.sfd_range[1] * 1e-3
        inrange = (Dvals >= Dkm_lo) & (Dvals <= Dkm_hi)
        lo = Dvals < Dkm_lo
        hi = Dvals > Dkm_hi
        t = 1.0
        Nvals = production.function(diameter=Dvals*1e3,age=t*1e3)
        Nvals *= 1e6 # convert from m^-2 to km^-2
        ax.plot(Dvals[inrange], Nvals[inrange], '-', color='black', linewidth=1.0, zorder=50)
        ax.plot(Dvals[lo], Nvals[lo], '-.', color='orange', linewidth=2.0, zorder=50)
        ax.plot(Dvals[hi], Nvals[hi], '-.', color='orange', linewidth=2.0, zorder=50)

        plt.tick_params(axis='y', which='minor')
        plt.tight_layout()
        plt.show()    
   
   
    def plot_sampled_csfd():
        from cratermaker.components.production.powerlaw import PowerLawProduction
        fig = plt.figure(1, figsize=(8, 7))
        ax = {'Power Law': fig.add_subplot(121),
            'NPF (Moon)': fig.add_subplot(122)}

        production = {
                    'Power Law': PowerLawProduction(impact_velocity_model = "Moon_MBA"),
                    'NPF (Moon)': NeukumProduction(version="Moon")
                    }
        
        target = Target("Moon")
        area = 4*np.pi*target.radius**2 
        age = 4100.0
        x_min = 1e0
        x_max = 1e5
        y_min = 1e0
        y_max = 1e4
        diameter_range = (2e3,10000e3) # Range of diameters to generate in m
        nD = 1000
        Dvals = np.logspace(np.log10(x_min), np.log10(x_max), num=nD)
        Nevaluations = 100
        for key in ax:
            ax[key].title.set_text(key)
            ax[key].set_xscale('log')
            ax[key].set_yscale('log')
            ax[key].set_ylabel('$\\mathregular{N_{>D}}$')
            ax[key].set_xlabel('Diameter (km)')
            ax[key].set_xlim(x_min, x_max)
            ax[key].set_ylim(y_min, y_max)
            ax[key].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
            ax[key].yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
            ax[key].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
            ax[key].xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
            # Plot the sampled values
            for i in range(Nevaluations):
                Dsampled, _ = production[key].sample(age=age, diameter_range=diameter_range, area=area, return_age=False)
                Dsampled = np.sort(Dsampled)[::-1]
                Nsampled = range(1,len(Dsampled)+1) 
                ax[key].plot(Dsampled*1e-3, Nsampled, '-', color='cornflowerblue', linewidth=2.0, zorder=50, alpha=0.2)
               
            # Plot the production function 
            Nvals = production[key].function(diameter=Dvals*1e3,age=age)
            Nvals *= area # convert from per unit area to total number
            ax[key].plot(Dvals, Nvals, '-', color='black', linewidth=3.0, zorder=50)

        plt.tick_params(axis='y', which='minor')
        plt.tight_layout()
        plt.show()
    
            
    # plot_npf_csfd()
    # plot_npf_N1_vs_T()
    # plot_npf_fit()    
    # plot_npf_proj_csfd()
    plot_sampled_csfd()

