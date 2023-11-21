import cratermaker
import numpy as np
import matplotlib.pyplot as plt
from cratermaker import craterscaling as cs    

sim = cratermaker.Simulation()
cs.transient_to_final(100e3,sim.target,sim.rng)

sim.add_projectile(diameter=5.0e3,velocity=25e3)
sim.projectile

sim.crater

sim.add_crater(diameter=79.0)
sim.crater

# Sample data
Dtrans_orig = np.random.uniform(5e3,100e3,size=1000)
Dfinal_new = []
Dtrans_new = []
for D in Dtrans_orig:
    sim.add_crater(transient_diameter=D)
    if sim.crater.transient_diameter > sim.crater.diameter:
        print("Something went wrong: Dt->Dt")
        print(f"Df/Dt: {sim.crater.diameter/sim.crater.transient_diameter}")    
    Dfinal = sim.crater.diameter
    Dfinal_new.append(Dfinal)
    sim.add_crater(diameter=Dfinal)
    Dtrans_new.append(sim.crater.transient_diameter)
    if sim.crater.transient_diameter > sim.crater.diameter:
        print("Something went wrong: Df->Dt")
        print(f"Df/Dt: {sim.crater.diameter/sim.crater.transient_diameter}")
   

Dfinal_orig = np.array(Dfinal_new)
Dtrans_new = np.array(Dtrans_new)