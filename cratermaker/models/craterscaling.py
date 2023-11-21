import numpy as np
from numpy.random import Generator
from ..core.target import Target
from ..core.projectile import Projectile
from ..core import montecarlo as mc
from scipy.optimize import root_scalar

def get_simple_to_complex_transition_factors(target: Target, rng: Generator=None):
    if not isinstance(target, Target):
        raise TypeError("target must be an instance of Target")
    if rng and not isinstance(rng, Generator):
        raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
    if rng is None:
        rng = np.random.default_rng() 
   
    # These terms are used to compute the ratio of the transient crater to simple crater size       
    simple_enlargement_mean = 0.84 # See Melosh (1989) pg. 129 just following eq. 8.2.1
    simple_enlargement_std = 0.04 # Just using Pike (1980) fig. 9 the crater depth varies by about the same amount on either side of the transition so this is a reasonable assumption
    
    # These terms are used in the exponent in the final rim radius/ simple crater radius vs  final radius / transition radius relationship
    # See Holsapple (1993) eq. 28
    final_exp_mean = 0.079    
    final_exp_std = 0.0001 # We add noise because this is nature and nature messy
    complex_enlargement_factor = 1.02
   
    # These terms are used to compute the transition diameter as a function of gravity
    # The transition values come from CTEM and are a synthesis of Pike (1980), Croft (1985), Schenk et al. (2004).
    if target.transition_scale_type == "silicate":
        simple_complex_exp = -1.0303 
        simple_complex_mean = 2*16533.8 
        simple_complex_std = 0.04
    elif target.transition_scale_type == "ice":
        simple_complex_exp = -1.22486
        simple_complex_mean = 2*3081.39
        simple_complex_std = 0.04
    
    # The nominal value will be used for determining the range of the "transitional" morphology type
    transition_nominal= simple_complex_mean * target.gravity**simple_complex_exp
    
    # Draw from a truncated normal distribution for each component of the model
    simple_enlargement_factor = 1.0 / mc.bounded_norm(simple_enlargement_mean, simple_enlargement_std)
    final_exp = mc.bounded_norm(final_exp_mean, final_exp_std)
    simple_complex_fac = simple_complex_mean * np.exp(rng.normal(loc=0.0,scale=simple_complex_std))
    transition_diameter = simple_complex_fac * target.gravity**simple_complex_exp
    factors = {
        "transition_diameter": transition_diameter,
        "transition_nominal" : transition_nominal,
        "simple_enlargement_factor" : simple_enlargement_factor,
        "complex_enlargement_factor" : complex_enlargement_factor,
        "final_exp" : final_exp,
        }
    return factors

def get_morphology_type(diameter, target: Target, rng: Generator=None, factors=None):
    if not isinstance(target, Target):
        raise TypeError("target must be an instance of Target")
    if rng and not isinstance(rng, Generator):
        raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
    if rng is None:
        rng = np.random.default_rng()    
   
    if factors is None: 
        f = get_simple_to_complex_transition_factors(target,rng)
    else:
        f = factors
        
    # Use the 1/2x to 2x the nominal value of the simple->complex transition diameter to get the range of the "transitional" morphology type. This is supported by: Schenk et al. (2004) and Pike (1980) in particular  
    transition_range = (0.5*f['transition_nominal'],2*f['transition_nominal'])
    
    if diameter < transition_range[0]:
        morphology_type = "simple" 
    elif diameter > transition_range[1]:
        morphology_type = "complex"
    else:
        # We'll uses the distance from the nominal transition diameter to set a probability of being either simple, complex, or transitional.
        if diameter < f['transition_nominal']:
            p = (f['transition_nominal'] - diameter)/(f['transition_nominal'] - transition_range[0])
            categories = ["simple","transitional"]
            prob = [p, 1.0-p] 
            morphology_type = rng.choice(categories,p=prob)
        else:
            p = (diameter - f['transition_nominal'])/(transition_range[1] - f['transition_nominal'])
            categories = ["complex","transitional"]
            prob = [p, 1.0-p] 
            morphology_type = rng.choice(categories,p=prob)                
    return morphology_type

def f2t_simple(Df,f):
    Dt = Df / f['simple_enlargement_factor']
    return Dt
    
def f2t_complex(Df,f):
    Dt = Df / (f['simple_enlargement_factor'] * f['complex_enlargement_factor']) * (Df / f['transition_diameter'])**-f['final_exp']    
    return Dt 

def final_to_transient(final_diameter, target: Target, rng: Generator=None):
    if not isinstance(target, Target):
        raise TypeError("target must be an instance of Target")
    if rng and not isinstance(rng, Generator):
        raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
    if rng is None:
        rng = np.random.default_rng() 
        
    f = get_simple_to_complex_transition_factors(target,rng)
    morphology_type = get_morphology_type(final_diameter, target, rng, f) 
    
    if morphology_type == "simple": 
        transient_diameter =  f2t_simple(final_diameter,f)
    else:
        transient_diameter = f2t_complex(final_diameter,f)

    return np.float64(transient_diameter), morphology_type


def transient_to_final(transient_diameter, target: Target, rng: Generator=None):
    if not isinstance(target, Target):
        raise TypeError("target must be an instance of Target")
    if rng and not isinstance(rng, Generator):
        raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
    if rng is None:
        rng = np.random.default_rng() 
    
    f = get_simple_to_complex_transition_factors(target,rng)
    # In the transition region, a particular transient crater diameter could be associate with simple, complex, or transitional crater morphologies.
    # Therefore we need to monte carlo our way into a solution to avoid biasing in favor of one or another in the transient->final computation
    
    # Invert the final -> transient functions for  each crater type
    final_diameter_simple = transient_diameter * f['simple_enlargement_factor']
    sol = root_scalar(lambda x, *args: f2t_complex(x, *args)-transient_diameter,bracket=(0.1*final_diameter_simple,10*final_diameter_simple), args=(f))
    final_diameter_complex = sol.root
    
    # Evaluate the potential morphology that this transient crater could be consistent with. If both potential diameter values are unambigusously is unambiguosuly simple or complex, go with that.
    # If there is disagreement, then we'll draw the answer from a hat and just check to make sure that final_diameter > transient_diameter 
    morphology_options = [get_morphology_type(final_diameter_simple,target,rng,f),get_morphology_type(final_diameter_complex,target,rng,f)]
    
    if len(set(morphology_options)) == 1: # We have agreement!
        morphology_type = morphology_options[0]
        if morphology_type == "simple":
            final_diameter = final_diameter_simple
        else:
            final_diameter = final_diameter_complex # this includes transitional types as well
    else: 
        if "simple" in morphology_options: # The disagreement is between simple/complex or simple/transitional
            if morphology_options[0] == "simple":
                sind = 0
                cind = 1 
            else:
                sind = 1
                cind = 0
                
            # Randomly draw a morphology based on weighting by whichever option is closest to the transition 
            is_simple = rng.random() < np.abs(final_diameter_complex - f['transition_diameter']) / np.abs(final_diameter_simple - final_diameter_complex)
            if is_simple:
                final_diameter = final_diameter_simple
                morphology_type = morphology_options[sind] 
            else:
                final_diameter = final_diameter_complex
                morphology_type = morphology_options[cind]
        else:
            final_diameter = final_diameter_complex
            morphology_type = rng.choice(morphology_options)
    
    
    return np.float64(final_diameter), morphology_type


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