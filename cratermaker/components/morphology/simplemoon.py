import os
import numpy as np
from numpy.random import Generator
import json
import math
from scipy import fft
from scipy.optimize import fsolve
from numpy.typing import NDArray, ArrayLike
from typing import Any
from cratermaker.core.target import Target
from cratermaker.utils.custom_types import FloatLike
from cratermaker.core.surface import Surface
from cratermaker.components.morphology import register_morphology_model, MorphologyModel
from cratermaker.utils.general_utils import parameter
from cratermaker import crater, ejecta

@register_morphology_model("simplemoon")
class SimpleMoon(MorphologyModel):
    """
    An operations class for computing the morphology of a crater based on its size and target properties.

    This class encapsulates the logic for altering the topography of the surface based on the crater properties.

    Parameters
    ----------
    crater : Crater
        The crater to be created.
    target : Target
        The target body for the impact simulation.
    ejecta_truncation : float, optional
        The relative distance from the rim of the crater to truncate the ejecta blanket, default is None, which will compute a 
        truncation distance based on where the ejecta thickness reaches a small value.         
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used.
    dorays : bool, optional
        A flag to determine if the ray pattern should be used instead of the homogeneous ejecta blanket, default is True.
    **kwargs : Any
        Additional keyword arguments to be passed to internal functions.
    """
    
    def __init__(self, 
                 crater,  
                 target: Target | None = None, 
                 ejecta_truncation: FloatLike | None = None,
                 rng: Generator | None = None,
                 dorays: bool = True,
                 **kwargs: Any 
                 ):
        
        self.crater = crater 
        self.target = target
        self.rng = rng
        self.ejecta_truncation = ejecta_truncation
        self.dorays = dorays

        # Set the morphology based on crater type
        diameter_m = self.crater.final_diameter
        diameter_km = diameter_m * 1e-3  # Convert to km for some models

        if self.crater.morphology_type in ["simple", "transitional"]:
            # A hybrid model between Pike (1977) and Fassett & Thomson (2014)
            self.rimheight = 0.043 * diameter_km**1.014 * 1e3  # Closer to Fassett & Thomson
            self.rimwidth = 0.257 * diameter_km**1.011 * 1e3   # Pike model
            self.floordepth = 0.224 * diameter_km**1.010 * 1e3 # Closer to Fassett & Thomson
            self.floor_diameter = 0.200 * diameter_km**1.143 * 1e3  # Fassett & Thomson for D~1km, Pike for D~20km

        elif self.crater.morphology_type in ["complex", "peakring", "multiring"]:
            # Following Pike (1977)
            self.rimheight = 0.236 * diameter_km**0.399 * 1e3  # Pike model
            self.rimwidth = 0.467 * diameter_km**0.836 * 1e3   # Pike model
            self.floordepth = 1.044 * diameter_km**0.301 * 1e3 # Pike model
            # Fassett & Thomson for D~1km, Pike for D~20km, but limited to 90% of diameter
            self.floor_diameter = min(0.187 * diameter_km**1.249 * 1e3, 0.9 * diameter_m)
            self.peakheight = 0.032 * diameter_km**0.900 * 1e3  # Pike model
        else:
            raise ValueError(f"Unknown morphology type: {self.morphology_type}")
            
        self.ejrim = 0.14 * (diameter_m * 0.5)**(0.74) # McGetchin et al. (1973) Thickness of ejecta at rim

       
    def __repr__(self):
        return (f"Morphology(morphology_type={self.crater.morphology_type}, diameter={self.crater.final_diameter}, "
                f"rimheight: {self.rimheight}, rimwidth: {self.rimwidth}, floordepth: {self.floordepth}, floor_diameter: {self.floor_diameter})") 
    

    def crater_profile(self, r: ArrayLike, r_ref: ArrayLike) -> np.float64:
        elevation = crater.profile(r,
                                   r_ref, 
                                   self.crater.final_diameter, 
                                   self.floordepth, 
                                   self.floor_diameter, 
                                   self.rimheight, 
                                   self.ejrim
                                )
        
        return np.array(elevation, dtype=np.float64)
    

    def ejecta_profile(self, r: ArrayLike) -> np.float64:
        elevation = ejecta.profile(r,
                                   self.crater.final_diameter, 
                                   self.ejrim
                                )
        elevation = np.array(elevation, dtype=np.float64)
        return elevation
   
    
    def ejecta_distribution(self, r: ArrayLike, theta: ArrayLike) -> np.float64:
        thickness = ejecta.distribution(r, theta,
                                       self.crater.final_diameter, 
                                       self.ejrim, 
                                       self.ejecta_truncation,
                                       self.dorays
                                    )
        thickness = np.array(thickness, dtype=np.float64)
        return thickness


    def ray_intensity(self, r: ArrayLike, theta: ArrayLike) -> np.float64:
        intensity = ejecta.ray_intensity(r, theta,
                                       self.crater.final_diameter, 
                                       self.ejrim, 
                                       self.ejecta_truncation,
                                    )
        intensity = np.array(intensity, dtype=np.float64)
        return intensity


    def compute_rmax(self, 
                     minimum_thickness: FloatLike,
                     feature: str = "ejecta") -> float:
        """
        Compute the maximum extent of the crater based on the minimum thickness of a feature, or the ejecta_truncation factor,
        whichever is smaller.
        
        Parameters
        ----------
        minimum_thickness : FloatLike
            The minimum thickness of the feature blanket in meters.
        feature : str, optional, default = "ejecta"
            The feature to compute the maximum extent. Either "crater" or "ejecta". If "crater" is chosen, the rmax is based
            on where the raised rim is smaller than minimum thickness. 
        Returns
        -------
        float
            The maximum extent of the crater or ejecta blanket in meters.
        """ 
        
        # Compute the reference surface for the crater 

        if feature == "ejecta":
            def _profile_invert(r):
                return self.ejecta_profile(r) - minimum_thickness
        elif feature == "crater":
            def _profile_invert(r):
                return self.crater_profile(r, np.zeros(1)) - minimum_thickness
        else:
            raise ValueError("Unknown feature type. Choose either 'crater' or 'ejecta'")
    
        # Get the maximum extent 
        rmax = fsolve(_profile_invert, x0=self.crater.final_radius*1.01)[0]
        
        if self.ejecta_truncation:
            rmax = min(rmax, self.ejecta_truncation * self.crater.final_radius)
            
        if feature == "crater":
            rmax = max(rmax, self.crater.final_radius)

        return float(rmax)


    def form_crater(self, 
                    surf: Surface,
                    **kwargs) -> None:
        """
        This method forms the interior of the crater by altering the elevation variable of the surface mesh.
        
        Parameters
        ----------
        surf : Surface
            The surface to be altered.
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions (not used here).
        """
        
        # Test if the crater is big enough to modify the surface
        rmax = self.compute_rmax(minimum_thickness=surf.smallest_length)
        region_surf = surf.extract_region(self.crater.location, rmax)
        if not region_surf: # The crater is too small to change the surface
            return
        crater_area = np.pi * rmax**2
        
        # Check to make sure that the face at the crater location is not smaller than the crater area
        if surf['face_areas'].isel(n_face=self.crater.face_index) > crater_area:
            return
        
        region_surf['node_crater_distance'], region_surf['face_crater_distance'] = region_surf.get_distance(self.crater.location)
        region_surf.get_reference_surface(self.crater.location, self.crater.final_radius)
        
        try:
            node_elevation = self.crater_profile(region_surf['node_crater_distance'].values, 
                                                 region_surf['reference_node_elevation'].values)
            surf['node_elevation'].loc[{'n_node': region_surf.uxgrid._ds["subgrid_node_indices"]}] = node_elevation
            
            face_elevation = self.crater_profile(region_surf['face_crater_distance'].values, 
                                                 region_surf['reference_face_elevation'].values)
            surf['face_elevation'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] = face_elevation
        except:
            print(self)
            raise ValueError("Something went wrong with this crater!")
                 
        return  


    def form_ejecta(self,
                    surf: Surface,
                    **kwargs) -> None:
        """
        This method forms the ejecta blanket around the crater by altering the elevation variable of the surface mesh.
       
        Parameters
        ----------
        surf : Surface
            The surface to be altered.
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions (not used here). 
        """
                 
        # Test if the ejecta is big enough to modify the surface
        rmax = self.compute_rmax(minimum_thickness=surf.smallest_length) 
        if not self.ejecta_truncation:
            self.ejecta_truncation = rmax / self.crater.final_radius
        region_surf = surf.extract_region(self.crater.location, rmax)
        if not region_surf: # The crater is too small to change the surface
            return
        ejecta_area = np.pi * rmax**2
        
        # Check to make sure that the face at the crater location is not smaller than the ejecta blanket area
        if surf['face_areas'].isel(n_face=self.crater.face_index) > ejecta_area:
            return                  
        
        region_surf['node_crater_distance'], region_surf['face_crater_distance'] = region_surf.get_distance(self.crater.location)
        region_surf['node_crater_bearing'], region_surf['face_crater_bearing']  = region_surf.get_initial_bearing(self.crater.location)
        
        try:
            node_thickness = self.ejecta_distribution(region_surf['node_crater_distance'].values, 
                                                     region_surf['node_crater_bearing'].values)
            surf['node_elevation'].loc[{'n_node': region_surf.uxgrid._ds["subgrid_node_indices"]}] += node_thickness
            
            face_thickness = self.ejecta_distribution(region_surf['face_crater_distance'].values, 
                                                     region_surf['face_crater_bearing'].values)
            surf['face_elevation'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] += face_thickness
            surf['ejecta_thickness'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] += face_thickness
            
            if self.dorays: 
                face_intensity = self.ray_intensity(region_surf['face_crater_distance'].values, 
                                                     region_surf['face_crater_bearing'].values)
                surf['ray_intensity'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] += face_intensity
        except:
            print(self)
            raise ValueError("Something went wrong with this crater!")
                 
        return  
    

    def form_secondaries(self,
                         surf: Surface,
                         **kwargs) -> None:
        """
        This method forms secondary craters around the primary crater. Currently it only generates the ray intensity function.
        Maximum ray length formula is from [1]_.
       
        Parameters
        ----------
        surf : Surface
            The surface to be altered.
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions (not used here). 
            
        References
        ----------
        .. [1] Elliott, J.R., Huang, Y.-H., Minton, D.A., Freed, A.M., 2018. The length of lunar crater rays explained using secondary crater scaling. Icarus 312, 231-246. https://doi.org/10.1016/j.icarus.2018.04.015

        """
        
        # Elliott et al. (2018) eq. 2
        A = 6.59 + self.rng.normal(loc=0.0,scale=1.45) 
        p = 1.27 + self.rng.normal(loc=0.0,scale=0.06)
        rmax = A * (self.crater.final_radius / 1000) **p * 1000
        if not self.ejecta_truncation:
            self.ejecta_truncation = rmax / self.crater.final_radius
        region_surf = surf.extract_region(self.crater.location, rmax)
        if not region_surf: # The crater is too small to change the surface
            return
        ray_area = np.pi * rmax**2
        
        # Check to make sure that the face at the crater location is not smaller than the ray effectt area
        if surf['face_areas'].isel(n_face=self.crater.face_index) > ray_area:
            return                  
        
        region_surf['node_crater_distance'], region_surf['face_crater_distance'] = region_surf.get_distance(self.crater.location)
        region_surf['node_crater_bearing'], region_surf['face_crater_bearing']  = region_surf.get_initial_bearing(self.crater.location)
        
        try:
            face_intensity = self.ray_intensity(region_surf['face_crater_distance'].values, 
                                                region_surf['face_crater_bearing'].values)
            surf['ray_intensity'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] += face_intensity
        except:
            print(self)
            raise ValueError("Something went wrong with this crater!")
                 
        return      


    def get_1D_power_spectral_density(self,feature,psd_coef,num_psd_component_effec=100) -> NDArray:
        """
        This method construct a 1D power spectral density.
        Coeffcients are from [1]_.
       
        Parameters
        ----------
        feature : string
            For a 1D feature, choose from ejecta, rim, and floor
        psd_coef : dict (.json)
            Coeffcients used to constract a 1D power spectral density
        num_psd_component_effec : int
            Only reconstrut the sine waves with wavelengths smaller than 2pi/num_psd_component_effec to improve computational efficiency

        References
        ----------
        .. [1] Du, J., Minton, D. A., Blevins, A. M., Fassett, C. I., & Huang, Y. H. (2024). Spectral analysis of the morphology of fresh lunar craters I: Rim crest, floor, and rim flank outlines. Journal of Geophysical Research: Planets, 129(11), e2024JE008357. https://doi.org/10.1029/2024JE008357

        """
        # read the psd_coef outside this function. this is temporary until we find a better way to read the psd_coef
        with open(os.path.join(os.pardir,"components","psd_coef.json")) as f:
            psd_coef = json.load(f)
        # ------------------------------------------------------------------------------------------------------------------
        num_psd_component=5000
        # ------------------------------------------------------------------------------------------------------------------
        if self.crater.final_diameter<psd_coef["1D"][feature]["slope_12"]["D_tie"]:
            slope_12 = psd_coef["1D"][feature]["slope_12"]["k1"] * self.crater.final_diameter + psd_coef["1D"][feature]["slope_12"]["b1"]
            slope_12_sigma = psd_coef["1D"][feature]["slope_12"]["sigma_1"]
        else:
            slope_12 = psd_coef["1D"][feature]["slope_12"]["k2"] * self.crater.final_diameter + psd_coef["1D"][feature]["slope_12"]["b2"]
            slope_12_sigma = psd_coef["1D"][feature]["slope_12"]["sigma_2"]
        if self.crater.final_diameter<psd_coef["1D"][feature]["bp2_x"]["D_tie"]:
            bp2_x = psd_coef["1D"][feature]["bp2_x"]["k1"] * self.crater.final_diameter + psd_coef["1D"][feature]["bp2_x"]["b1"]
            bp2_x_sigma = psd_coef["1D"][feature]["bp2_x"]["sigma_1"]
        else:
            bp2_x = psd_coef["1D"][feature]["bp2_x"]["k2"] * self.crater.final_diameter + psd_coef["1D"][feature]["bp2_x"]["b2"]
            bp2_x_sigma = psd_coef["1D"][feature]["bp2_x"]["sigma_2"]
        if self.crater.final_diameter<psd_coef["1D"][feature]["bp2_y"]["D_tie"]:
            bp2_y = psd_coef["1D"][feature]["bp2_y"]["k1"] * self.crater.final_diameter + psd_coef["1D"][feature]["bp2_y"]["b1"]
            bp2_y_sigma = psd_coef["1D"][feature]["bp2_y"]["sigma_1"]
        else:
            bp2_y = psd_coef["1D"][feature]["bp2_y"]["k2"] * self.crater.final_diameter + psd_coef["1D"][feature]["bp2_y"]["b2"]
            bp2_y_sigma = psd_coef["1D"][feature]["bp2_y"]["sigma_2"]
        if self.crater.final_diameter<psd_coef["1D"][feature]["bp3_y"]["D_tie"]:
            bp3_y = psd_coef["1D"][feature]["bp3_y"]["k1"] * self.crater.final_diameter + psd_coef["1D"][feature]["bp3_y"]["b1"]
            bp3_y_sigma = psd_coef["1D"][feature]["bp3_y"]["sigma_1"]
        else:
            bp3_y = psd_coef["1D"][feature]["bp3_y"]["k2"] * self.crater.final_diameter + psd_coef["1D"][feature]["bp3_y"]["b2"]
            bp3_y_sigma = psd_coef["1D"][feature]["bp3_y"]["sigma_2"]
        if self.crater.final_diameter<psd_coef["1D"][feature]["bp4_y"]["D_tie"]:
            bp4_y = psd_coef["1D"][feature]["bp4_y"]["k1"] * self.crater.final_diameter + psd_coef["1D"][feature]["bp4_y"]["b1"]
            bp4_y_sigma = psd_coef["1D"][feature]["bp4_y"]["sigma_1"]
        else:
            bp4_y = psd_coef["1D"][feature]["bp4_y"]["k2"] * self.crater.final_diameter + psd_coef["1D"][feature]["bp4_y"]["b2"]
            bp4_y_sigma = psd_coef["1D"][feature]["bp4_y"]["sigma_2"]
        slope_12 = np.random.normal(slope_12, slope_12_sigma)
        bp2_x = np.random.normal(bp2_x, bp2_x_sigma)
        bp2_y = np.random.normal(bp2_y, bp2_y_sigma)
        bp3_y = np.random.normal(bp3_y, bp3_y_sigma)
        bp4_y = np.random.normal(bp4_y, bp4_y_sigma)
        psd_sigma=psd_coef["1D"][feature]["psd_sigma"]
        # ------------------------------------------------------------------------------------------------------------------
        bp4_x = math.log10(2*math.pi)
        bp3_x =  math.log10(10**bp4_x / 2)
        interval=2*math.pi/num_psd_component
        freq = fft.fftfreq(num_psd_component, interval)
        wavelength = 1 / freq[1:num_psd_component//2]
        psd = np.zeros((num_psd_component//2 - 1, 2))
        psd[:, 0] = wavelength
        # ----------------------------------------------------------------------------------------------------------------------
        bp2_x_index = int(  10**bp4_x/10**bp2_x-1)
        k_23=(bp3_y-bp2_y)/(bp3_x-bp2_x)
        b_23 = bp3_y-k_23*bp3_x
        k_12=slope_12
        b_12 = bp2_y-k_12*bp2_x
        #-----------------------------------------------------------------------------------------------------------------------
        psd[0,1] =10**bp4_x
        psd[1,1]=10**bp3_y
        psd[2:bp2_x_index+1,1] =10**(k_23*np.log10(psd[2:bp2_x_index+1, 0])+b_23)
        psd[bp2_x_index+1:, 1] =10**(k_12*np.log10(psd[bp2_x_index+1:, 0])+b_12)
        # ----------------------------------------------------------------------------------------------------------------------
        psd=psd[:num_psd_component_effec]
        psd_log = np.log10(psd[:, 1])
        psd_log += np.random.normal(0, psd_sigma, psd_log.shape) 
        psd[:, 1] = 10 ** psd_log
        psd=np.flipud(psd)
        return psd


    def get_2D_power_spectral_density(self,feature,psd_coef,ejecta_radius_norm,floor_radius_norm,max_effec_freq) -> tuple[NDArray,NDArray,NDArray,NDArray,NDArray,NDArray]:
        """
        This method constructs a 2D power spectral density.
        Coeffcients are from [1]_.
       
        Parameters
        ----------
        feature : string
            For a 2D feature, choose from ejecta, wall, and floor
        psd_coef : dict (.json)
            Coeffcients used to constract a 2D power spectral density
        ejecta_radius_norm : float
            The radius of the continuous ejecta normalized by the crater radius
        floor_radius_norm : float
            The radius of the floor normalized by the crater radius
        max_effec_freq : int
            Only reconstrut the sine waves with frequencies smaller than a maxmium effective frequency to improve computational efficiency 

        References
        ----------
        .. [1] Du, J., Minton, D. A., Blevins, A. M., Fassett, C. I., & Huang, Y. H. (2025). Spectral Analysis of the Morphology of Fresh Lunar Craters II: Two-Dimensional Surface Elevations of the Continuous Ejecta, Wall, and Floor. Journal of Geophysical Research: Planets

        """
         # read the psd_coef outside this function. this is temporary until we find a better way to read the psd_coef
        with open(os.path.join(os.pardir,"components","psd_coef.json")) as f:
            psd_coef = json.load(f) # read in from _init_?
        # ------------------------------------------------------------------------------------------------------------------
        if feature == "ejecta" : max_effec_freq = 40
        if feature == "wall" : max_effec_freq = 40
        if feature == "floor" : max_effec_freq = 100
        # ------------------------------------------------------------------------------------------------------------------
        if feature=="ejecta":
            theta_min = 0
            theta_max = 2 * math.pi
            r_min=1
            r_max=ejecta_radius_norm
            area=(r_max - r_min) * (theta_max - theta_min)
            theta_number = 2000
            r_number = 400
        if feature == "wall":
            theta_min = 0
            theta_max = 2 * math.pi
            r_min = floor_radius_norm
            r_max = 1
            area= (r_max - r_min) * (theta_max - theta_min)
            theta_number = 1000
            r_number = 200
        if feature == "floor":
            r_min = -floor_radius_norm/math.sqrt(2)
            r_max = floor_radius_norm/math.sqrt(2)
            theta_min = r_min
            theta_max = r_max
            area=(r_max - r_min) ** 2
            theta_number = 512
            r_number = 512
        # ------------------------------------------------------------------------------------------------------------------
        side_r=np.linspace(0, r_max-r_min, num=r_number, endpoint=False)
        side_theta=np.linspace(0, theta_max-theta_min, num=theta_number, endpoint=False)
        grid_theta, grid_r = np.meshgrid( side_theta,side_r)
        side_freq_theta = np.fft.fftfreq(theta_number, side_theta[1] - side_theta[0])
        side_freq_r = np.fft.fftfreq(r_number, side_r[1] - side_r[0])
        freq_theta, freq_r = np.meshgrid(np.fft.fftshift(side_freq_theta), np.fft.fftshift(side_freq_r))
        freq_theta_quadrant = freq_theta [r_number // 2:, theta_number // 2:]
        freq_r_quadrant = freq_r  [r_number // 2:, theta_number // 2:]
        # ------------------------------------------------------------------------------------------------------------------
        if self.crater.final_diameter<psd_coef["2D"][feature]["p_max"]["D_tie"]:
            p_max = psd_coef["2D"][feature]["p_max"]["k1"] * self.crater.final_diameter + psd_coef["2D"][feature]["p_max"]["b1"]
            p_max_sigma = psd_coef["2D"][feature]["p_max"]["sigma_1"]
        else:
            p_max = psd_coef["2D"][feature]["p_max"]["k2"] * self.crater.final_diameter + psd_coef["2D"][feature]["p_max"]["b2"]
            p_max_sigma = psd_coef["2D"][feature]["p_max"]["sigma_2"]
        if self.crater.final_diameter<psd_coef["2D"][feature]["p_diff"]["D_tie"]:
            p_diff = psd_coef["2D"][feature]["p_diff"]["k1"] * self.crater.final_diameter + psd_coef["2D"][feature]["p_diff"]["b1"]
            p_diff_sigma = psd_coef["2D"][feature]["p_diff"]["sigma_1"]
        else:
            p_diff = psd_coef["2D"][feature]["p_diff"]["k2"] * self.crater.final_diameter + psd_coef["2D"][feature]["p_diff"]["b2"]
            p_diff_sigma = psd_coef["2D"][feature]["p_diff"]["sigma_2"]
        if self.crater.final_diameter<psd_coef["2D"][feature]["nu_fall"]["D_tie"]:
            nu_fall = psd_coef["2D"][feature]["nu_fall"]["k1"] * self.crater.final_diameter + psd_coef["2D"][feature]["nu_fall"]["b1"]
            nu_fall_sigma = psd_coef["2D"][feature]["nu_fall"]["sigma_1"]
        else:
            nu_fall = psd_coef["2D"][feature]["nu_fall"]["k2"] * self.crater.final_diameter + psd_coef["2D"][feature]["nu_fall"]["b2"]
            nu_fall_sigma = psd_coef["2D"][feature]["nu_fall"]["sigma_2"]
        if self.crater.final_diameter<psd_coef["2D"][feature]["psd_sigma"]["D_tie"]:
            psd_sigma = psd_coef["2D"][feature]["psd_sigma"]["k1"] * self.crater.final_diameter + psd_coef["2D"][feature]["psd_sigma"]["b1"]
            psd_sigma_sigma = psd_coef["2D"][feature]["psd_sigma"]["sigma_1"]
        else:
            psd_sigma = psd_coef["2D"][feature]["psd_sigma"]["k2"] * self.crater.final_diameter + psd_coef["2D"][feature]["psd_sigma"]["b2"]
            psd_sigma_sigma = psd_coef["2D"][feature]["psd_sigma"]["sigma_2"]
        p_max = np.random.normal(p_max, p_max_sigma)
        p_diff = np.random.normal(p_diff, p_diff_sigma)
        nu_fall = np.random.normal(nu_fall, nu_fall_sigma)
        psd_sigma = np.random.normal(psd_sigma, psd_sigma_sigma)
        # ------------------------------------------------------------------------------------------------------------------
        freq_r_matrix= np.sqrt(freq_theta_quadrant**2+freq_r_quadrant**2)
        freq_theta_matrix = np.arctan2(freq_r_quadrant,freq_theta_quadrant)
        if feature=="ejecta":
            if self.crater.final_diameter<psd_coef["2D"][feature]["E_rad"]["D_tie"]:
                E_rad = psd_coef["2D"][feature]["E_rad"]["k1"] * self.crater.final_diameter + psd_coef["2D"][feature]["E_rad"]["b1"]
                E_rad_sigma = psd_coef["2D"][feature]["E_rad"]["sigma_1"]
            else:
                E_rad = psd_coef["2D"][feature]["E_rad"]["k2"] * self.crater.final_diameter + psd_coef["2D"][feature]["E_rad"]["b2"]
                E_rad_sigma = psd_coef["2D"][feature]["E_rad"]["sigma_2"]
            E_rad = np.random.normal(E_rad, E_rad_sigma)
            psd_log_quadrant=p_diff/np.sqrt(1-(E_rad*np.sin(freq_theta_matrix))**2)*(np.exp(-np.sqrt(freq_r_matrix/nu_fall))-1)+p_max
        if feature=="wall":
            if self.crater.final_diameter<psd_coef["2D"][feature]["E_circ"]["D_tie"]:
                E_circ = psd_coef["2D"][feature]["E_circ"]["k1"] * self.crater.final_diameter + psd_coef["2D"][feature]["E_circ"]["b1"]
                E_circ_sigma = psd_coef["2D"][feature]["E_circ"]["sigma_1"]
            else:
                E_circ = psd_coef["2D"][feature]["E_circ"]["k2"] * self.crater.final_diameter + psd_coef["2D"][feature]["E_circ"]["b2"]
                E_circ_sigma = psd_coef["2D"][feature]["E_circ"]["sigma_2"]
            E_circ = np.random.normal(E_circ, E_circ_sigma)
            psd_log_quadrant=p_diff/np.sqrt(1-(E_circ*np.cos(freq_theta_matrix))**2)*(np.exp(-np.sqrt(freq_r_matrix/nu_fall))-1)+p_max
        if feature=="floor":
            psd_log_quadrant = p_diff  * (np.exp(-np.sqrt(freq_r_matrix / nu_fall)) - 1) + p_max
        psd_log_quadrant_1st=np.random.normal(psd_log_quadrant,psd_sigma)
        psd_log_quadrant_2nd = np.random.normal(psd_log_quadrant, psd_sigma)
        # ------------------------------------------------------------------------------------------------------------------
        psd_log = np.zeros((r_number, theta_number))
        psd_log[r_number//2:, theta_number // 2:] =  psd_log_quadrant_1st
        psd_log[r_number//2+1:, 1:theta_number // 2] = np.fliplr(psd_log_quadrant_2nd[1:,1:])
        psd_log[ 1:r_number // 2, theta_number // 2+1:] = np.flipud(psd_log_quadrant_2nd[1:, 1:])
        psd_log[1:r_number // 2, 1:theta_number // 2] = np.fliplr(np.flipud(psd_log_quadrant_1st[1:, 1:]))
        #---------------------------------------------------------------------------------------------------------------
        psd_log[1:r_number // 2, theta_number // 2] = np.flip(psd_log[r_number // 2+1:, theta_number // 2])
        psd_log[r_number // 2, 1:theta_number // 2] = np.flip(  psd_log[r_number // 2, theta_number // 2+1:] )
        #---------------------------------------------------------------------------------------------------------------
        psd_log[0,theta_number//2+1:]=psd_log[1,theta_number//2+1:]
        psd_log[r_number // 2 + 1:, 0] = psd_log[r_number // 2 + 1:, 1]
        #---------------------------------------------------------------------------------------------------------------
        psd_log[0, 1:theta_number // 2] = np.flip(psd_log[0,theta_number//2+1:])
        psd_log[1:r_number // 2:, 0] = np.flip(psd_log[r_number // 2 + 1:,0] )
        # --------------------------------------------------------------------------------------------------------------
        psd_log[0, 0]=psd_log[1, 1]
        psd_log[0, theta_number//2] =psd_log[1, theta_number//2]
        psd_log[r_number//2, 0] = psd_log[r_number//2, 1]
        # ------------------------------------------------------------------------------------------------------------------
        psd = 10 ** psd_log
        phase = np.random.uniform(-math.pi, math.pi, size=(r_number, theta_number))
        if feature=="ejecta" or feature=="wall" : psd[:,theta_number//2]=np.zeros((r_number))
        # ------------------------------------------------------------------------------------------------------------------
        freq_theta_delta=(freq_theta.max()-freq_theta.min())/theta_number
        freq_r_delta = (freq_r.max() - freq_r.min()) / r_number
        freq_theta_num=int(max_effec_freq/freq_theta_delta)
        freq_r_num = int(max_effec_freq / freq_r_delta)
        freq_theta = freq_theta[r_number//2- freq_r_num:r_number//2+ freq_r_num,theta_number//2-freq_theta_num:theta_number//2+ freq_theta_num]
        freq_r = freq_r[r_number//2- freq_r_num:r_number//2+ freq_r_num,theta_number//2-freq_theta_num:theta_number//2+ freq_theta_num]
        psd = psd[r_number//2- freq_r_num:r_number//2+ freq_r_num,theta_number//2-freq_theta_num:theta_number//2+ freq_theta_num]
        phase = phase[r_number//2- freq_r_num:r_number//2+ freq_r_num,theta_number//2-freq_theta_num:theta_number//2+ freq_theta_num]
        # ------------------------------------------------------------------------------------------------------------------
        psd=np.sqrt(psd*area)
        return psd,phase,freq_theta,freq_r,grid_theta,grid_r  

    
    @property
    def rimheight(self) -> float:
        """
        The height of the crater rim in meters.
        
        Returns
        -------
        float
        """
        return self._rimheight
    
    @rimheight.setter
    def rimheight(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("rimheight must be of type FloatLike")
        self._rimheight = float(value)
        
    @property
    def rimwidth(self) -> float:
        """
        The width of the crater rim in meters.
        
        Returns
        -------
        float
        """
        return self._rimwidth
    
    @rimwidth.setter
    def rimwidth(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("rimwidth must be of type FloatLike") 
        self._rimwidth = float(value)
        
    @property
    def peakheight(self) -> float:
        """
        The height of the central peak in meters.
        
        Returns
        -------
        float
        """
        return self._peakheight
    
    @peakheight.setter
    def peakheight(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("peakheight must be of type FloatLike") 
        self._peakheight = float(value)

    @property
    def floor_diameter(self) -> float:
        """
        The diameter of the crater floor in meters.
        
        Returns
        -------
        float
        """
        return self._floor_diameter
    
    @floor_diameter.setter
    def floor_diameter(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("floor_diameter must be of type FloatLike")
        self._floor_diameter = float(value)
        
    @property
    def floordepth(self) -> float:
        """
        Return the depth of the crater floor in m
        
        Returns
        -------
        float
        """
        return self._floordepth
    
    @floordepth.setter
    def floordepth(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("floordepth must be of type FloatLike")
        self._floordepth = float(value)
        
    @property
    def ejrim(self) -> float:
        """
        The thickness of ejecta at the rim in m.
        
        Returns
        -------
        float
        """
        return self._ejrim
    
    @ejrim.setter
    def ejrim(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("ejrim must be of type FloatLike")
        self._ejrim = float(value)
        
    @property
    def crater(self):
        """
        The crater to be created.
        
        Returns
        -------
        Crater
        """ 
        return self._crater
    
    @crater.setter
    def crater(self, value):
        from cratermaker import Crater
        if value is not None and not isinstance(value, Crater):
            raise TypeError("crater must be an instance of Crater")
        self._crater = value
        return 
        
    @property
    def target(self):
        """
        The target body for the impact.
        
        Returns
        -------
        Target
        """ 
        return self._target
    
    @target.setter
    def target(self, value):
        if value is None:
            self._target = Target(name="Moon")
            return
        if not isinstance(value, Target):
            raise TypeError("target must be an instance of Target")
        self._target = value
        return        
        
    @property
    def rng(self):
        """
        A random number generator instance.
        
        Returns
        -------
        Generator
        """ 
        return self._rng
    
    @rng.setter
    def rng(self, value):
        if not isinstance(value, Generator) and value is not None:
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        self._rng = value or np.random.default_rng()   
       
    @parameter
    def ejecta_truncation(self) -> float:
        """
        The radius at which the crater is truncated relative to the crater radius.
        
        Returns
        -------
        float or None
        """
        return self._ejecta_truncation 
    
    @ejecta_truncation.setter
    def ejecta_truncation(self, value: FloatLike | None): 
        if value:
            if not isinstance(value, FloatLike):
                raise TypeError("truction_radius must be of type FloatLike")
            self._ejecta_truncation = float(value)
        else:
            self._ejecta_truncation = None
            
    @parameter
    def dorays(self) -> bool:
        """
        A flag to determine if the ray pattern should be used instead of the homogeneous ejecta blanket.
        
        Returns
        -------
        bool
        """
        return self._dorays
    
    @dorays.setter
    def dorays(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("dorays must be of type bool")
        self._dorays = value
