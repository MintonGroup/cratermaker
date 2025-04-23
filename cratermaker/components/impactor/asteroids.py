import numpy as np
from numpy.typing import NDArray
from typing import Any
from cratermaker.utils import montecarlo as mc
from cratermaker.components.impactor import ImpactorModel, register_impactor_model

@register_impactor_model("asteroids")
class AsteroidImpactors(ImpactorModel):
    """
    An operations class for computing the impactor properties of an asteroid source population.

    Parameters
    ----------
    **kwargs : Any
        Additional keyword arguments to be passed to internal functions.
    """
    
    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)



    @property
    def vertical_velocity(self):
        """Get the impact velocity model name."""
        if self._vertical_velocity is None:
            self.vertical_velocity = None
        return self._vertical_velocity


    @vertical_velocity.setter
    def vertical_velocity(self, value):
        """"
        Sets the vertical component of the impact velocity. If none is provided, it will look at the target name, and if it matches a known body, we will draw from a distribution using the predefined mean velocityvalue. Otherwise, it will raise an error.
            
        Notes
        ----- 
        Mean velocities for terrestrial planets and the Moon are based on analysis of simulations of main-belt derived asteroids from Minton & Malhotra (2010) [1]_  and Yue et al. (2013) [2]_. Mean velocities for the asteroids are from Bottke et al. (1994) [3]_.
        
        References
        ----------
        .. [1] Minton, D.A., Malhotra, R., 2010. Dynamical erosion of the asteroid belt and implications for large impacts in the inner Solar System. Icarus 207, 744-757. https://doi.org/10.1016/j.icarus.2009.12.008
        .. [2] Yue, Z., Johnson, B.C., Minton, D.A., Melosh, H.J., Di, K., Hu, W., Liu, Y., 2013. Projectile remnants in central peaks of lunar impact craters. Nature Geosci 6, 435 EP-. https://doi.org/10.1038/ngeo1828
        .. [3] Bottke, W.F., Nolan, M.C., Greenberg, R., Kolvoord, R.A., 1994. Velocity distributions among colliding asteroids. Icarus 107, 255-268. https://doi.org/10.1006/icar.1994.1021
        """

        if value is None: 
            pmv = self.mean_velocity
            vencounter_mean = np.sqrt(pmv**2 - self.target.escape_velocity**2)
            vencounter = mc.get_random_velocity(vencounter_mean, rng=self.rng)
            pv = np.sqrt(vencounter**2 + self.target.escape_velocity**2)
            pang = mc.get_random_impact_angle(rng=self.rng)
            self._vertical_velocity = pv * np.sin(np.deg2rad(pang))
        elif isinstance(value, (int, float)):
            if value < 0:
                raise ValueError("vertical_velocity must be a positive number")
            self._vertical_velocity = float(value)
        else: 
            raise TypeError("vertical_velocity must be a numeric value or None") 

        return

    @property
    def mean_velocity(self):
        """
        The mean velocity of the projectile in m/s.
        
        Returns
        -------
        float 
        """
        predefined_models = ['Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'MBA']
        predefined_velocities = [41100.0, 29100.0, 24600.0, 22100.0, 10700.0, 5300.0]
        predefined = dict(zip(predefined_models, predefined_velocities))
        if self.target.name in predefined_models:
            pmv = float(predefined[self.target.name])
        elif self.target.name in ["Ceres", "Vesta"]:
            pmv = float(predefined["MBA"])
        else:
            pmv = None
        return pmv
