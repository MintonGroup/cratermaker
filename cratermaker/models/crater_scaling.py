import numpy as np
from numpy.random import Generator
from ..core.target import Target
from ..core.projectile import Projectile

def get_simple_to_complex_transition_factors(target: Target, rng: Generator=None):
   if not isinstance(target, Target):
      raise TypeError("target must be an instance of Target")
   if rng and not isinstance(rng, Generator):
      raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
   if rng is None:
      rng = np.random.default_rng() 
        
   simple_enlargement_factor = 0.84 # See Melosh (1989) pg. 129 just following eq. 8.2.1
   transition_exp_mean = 0.15 # See Croft (1985) eq. 9
   final_exp_mean = 0.85      # See Croft (1985) eq. 9
   # Standard deviations for the exponents
   transition_exp_std = 0.04
   final_exp_std = 0.04
   
   # The transition values come from CTEM and are a synthesis of Pike (1980), Croft (1985), Schenk et al. (2004).
   if target.transition_scale_type == "silicate":
      simple_complex_exp = -1.0303 
      simple_complex_mean = 2*16533.8 
      simple_complex_std = 0.05
   elif target.transition_scale_type == "ice":
      simple_complex_exp = -1.22486
      simple_complex_mean = 2*3081.39
      simple_complex_std = 0.05
   

   # The nominal value will be used for determining the range of the "transitional" morphology type
   transition_nominal= simple_complex_mean * target.gravity**simple_complex_exp
   
   # Draw from a normal distribution for each exponent
   transition_exp = rng.normal(transition_exp_mean, transition_exp_std)
   final_exp = rng.normal(final_exp_mean, final_exp_std)
   simple_complex_fac = simple_complex_mean * np.exp(rng.normal(scale=simple_complex_std))
   transition_diameter = simple_complex_fac * target.gravity**simple_complex_exp
   return transition_diameter, transition_nominal, simple_enlargement_factor, transition_exp, final_exp


def final_to_transient(final_diameter, target: Target, rng: Generator=None):
   if not isinstance(target, Target):
      raise TypeError("target must be an instance of Target")
   if rng and not isinstance(rng, Generator):
      raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
   if rng is None:
      rng = np.random.default_rng()    
   transition_diameter, _, simple_enlargement_factor, transition_exp, final_exp = get_simple_to_complex_transition_factors(target,rng) 
   
   if final_diameter < transition_diameter:
      transient_diameter = simple_enlargement_factor * final_diameter  # Simple crater scaling
   else:
      transient_diameter = 1e3 * (transition_diameter * 1e-3)**transition_exp * (final_diameter * 1e-3)**final_exp # Complex crater scaling (in m)
   
   return np.float64(transient_diameter)


def transient_to_final(transient_diameter, target: Target, rng: Generator=None):
   if not isinstance(target, Target):
      raise TypeError("target must be an instance of Target")
   if rng and not isinstance(rng, Generator):
      raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
   if rng is None:
      rng = np.random.default_rng() 
           
   transition_diameter, _, simple_enlargement_factor, transition_exp, final_exp = get_simple_to_complex_transition_factors(target,rng)
   
   final_diameter_simple = transient_diameter / simple_enlargement_factor
   final_diameter_complex = 1e3 * ((1e-3 * transient_diameter) / (transition_diameter * 1e-3)**transition_exp)**(1.0 / final_exp)
   
   if final_diameter_simple < transition_diameter and final_diameter_complex < transition_diameter: # Unambiguosly simple
      final_diameter = final_diameter_simple
   elif final_diameter_simple > transition_diameter and final_diameter_complex > transition_diameter: # Unambiguously complex
      final_diameter = final_diameter_complex
   else: # Could be either complex or simple. We'll just draw which one from a hat weighted in the direction of whichever size is closest to the transition
      is_simple = rng.random() < np.abs(final_diameter_complex - transition_diameter) / np.abs(final_diameter_simple - final_diameter_complex)
      if is_simple:
         final_diameter = final_diameter_simple
      else:
         final_diameter = final_diameter_complex
      
   return np.float64(final_diameter)


def projectile_to_transient(projectile: Projectile, target: Target, rng: Generator=None):
   if not isinstance(projectile, Projectile):
      raise TypeError("target must be an instance of Target")
   if not isinstance(target, Target):
      raise TypeError("target must be an instance of Target")
   if rng and not isinstance(rng, Generator):
      raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
   if rng is None:
      rng = np.random.default_rng() 
      
   # Compute some auxiliary quantites
   projectile.mass = 4.0/3.0 * np.pi * projectile.density * (projectile.radius)**3
   mu = target.material.mu
   kv = target.material.K1
   c1 = 1.0 + 0.5 * mu
   c2 = (-3 * mu)/(2.0 + mu)

   # Find dimensionless quantities
   pitwo = (target.gravity * projectile.radius)/(projectile.vertical_velocity**2)
   pithree = target.material.Ybar / (target.material.density * (projectile.vertical_velocity**2))
   pifour = target.material.density / projectile.density
   pivol = kv * ((pitwo * (pifour**(-1.0/3.0))) + (pithree**c1))**c2
   pivolg = kv * (pitwo * (pifour**(-1.0/3.0)))**c2
   
   # find transient crater volume and radii (depth = 1/3 diameter)
   cvol = pivol * (projectile.mass / target.material.density)
   cvolg = pivolg * (projectile.mass / target.material.density)
   transient_radius = (3 * cvol / np.pi)**(1.0/3.0)
   transient_radius_gravscale = (3 * cvolg / np.pi)**(1.0/3.0)
    
   return transient_radius, transient_radius_gravscale