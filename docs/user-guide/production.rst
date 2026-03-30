.. currentmodule:: cratermaker.production

.. ipython:: python
    :okwarning:
    :suppress:

    from cratermaker import cleanup
    cleanup()

.. image:: ../_images/production_icon.svg
    :alt: Production
    :align: center
    :width: 300px
    :class: dark-light



.. _ug-production:

Production
==========

Cratermaker's production module provides a robust way to compute the production function for craters and projectiles. The production function is defined as the cumulative number of craters greater than a given diameter per unit m² surface area.  There are currently two types of production funcitons that are available in the |Production| class: "neukum" and "powerlaw". The default production model depends on the value of the target. If the target is one of the inner planets (or Vesta or Ceres), the default will be "neukum". For all other bodies, it is "powerlaw".  The "neukum" production function comes in three different versions: "moon", "mars", or "projectile."  The "projectile" version is used unless the target is either the Moon or Mars. The default Cratermaker target body is the Moon, and therefore the default production function is "neukum" with the "moon" model. We can create a production function object and inspect its values using the following code:

.. ipython:: python
   :okwarning:

   from cratermaker import Production
   production = Production.maker()
   print(production)

The output shows that by default, the production function is the Neukum model for the Moon, which is described in Neukum, Ivanov, and Hartmann (2001) [#]_. It also reports that the "Generator type" is "crater." A production model can either a crater generator or a projectile generator, which means that the diameter values that are returned from the production function represent either the crater or projectile. The other crater-based model is the Mars version, which is based on the work of Ivanov (2001) [#]_. The projectile model is based on the work of Ivanov, Neukum, and Wagner (2001) [#]_.

The Neukum Moon model has a valid range of ages over which it is valid, which is also reported in the output. The valid range of ages is between 0 and 4.5 Ga, though caution must be used when interpreting ages older than about 3.9 Ga, as these are poorly calibrated. The total crater number density of the most ancient terrains on the Moon correspond to a model age of 4.31 Ga, though the actual crust could potentially be must older than that. Nevertheless, the Neukum model is a well-established and widely used model for inner solar system crater chronology.

The Production class has two primary methods that are used to compute the production function. The first is |production.function|, which computes the expected cumulative number density of craters greater than a given diameter of a surface of a given age. 
The second is |production.sample|, which samples crater diameters and ages from the production function using Monte Carlo methods.

Production function
-------------------

|production.function| returns the cumulative size-frequency distribution (CSFD) of craters over a given age range and crater diameter. It takes the following arguments:

- **diameter**: Crater diameter(s) in units of meters to compute corresponding cumulative number density value.
- **time_start**: Starting time in units of My relative to the present, used to compute the CSFD. Default is ``1.0``, corresponding to 1 Ma. Alternatively `age` can be used as an alias for `time_start`.
- **time_end**: Ending time in units of My relative to the present, also used to compute the CSFD. Default is ``0.0``, which corresponds to the present day.

Example: Using Production.function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For this example, we are going to use |production.function| to find the cumulative number of craters greater than 1 km² to form on a surface in 3 Gy using a power law with a slope of -2. 

.. ipython:: python
   :okwarning:

   from cratermaker import Production 
   production = Production.maker("powerlaw",slope=-2.0) 
   diameter = 1000.0 
   age = 3000.0  
   n = production.function(diameter=diameter,age=age)  
   print(f"{n*1e6:.2f} D>{diameter*1e-3:.0f} km craters per km² in {age*1e-3:.0f} Gy")

Example: Plot a power law production function with a slope of -3 for craters 1 m to 100 km over 1 Ga
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. ipython:: python
   :okwarning:

   from cratermaker import Production
   import matplotlib.pyplot as plt
   import numpy as np

   fig, ax = plt.subplots()
   ax.set_xlabel('D (m)')
   ax.set_ylabel('$\\mathregular{N_{>D}}$')

   production = Production.maker("powerlaw", slope=-3.0)
   diameters = np.arange(1, 1e6)
   pf = production.function(diameter=diameters, age=1000.0)

   ax.loglog(diameters, pf)
   @savefig PowerLawExample.png width=4in
   plt.show()


Example: Using Production.age_from_D_N
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|production.age_from_D_N| returns the age in My for a given crater number density and diameter. It has the following parameters:

- **diameter**: Diameter of the crater in m
- **cumulative_number_density**: Number density of craters per m² surface area greater than the input diameter.

For this example, we are going to use |production.age_from_D_N| to plot the age in Ma for a number density of 10\ :sup:`-6` craters per m² from 1 m to 1 km in diameter. 


.. ipython:: python
   :okwarning:

   from cratermaker import Production
   import matplotlib.pyplot as plt
   import numpy as np

   production = Production.maker("powerlaw", slope=-2.5)
   diameters = np.logspace(0, 3)
   cumulative_number_density = np.full_like(diameters, 1e-6)
   age = production.age_from_D_N(diameter=diameters, cumulative_number_density=cumulative_number_density)

   def plot_inv(diameter, age):
      fig, ax = plt.subplots()
      ax.set_xlabel('D (m)')
      ax.set_ylabel('Age (Ma)')
      ax.loglog(diameters, age)
      ax.grid(True)
      return fig

   fig = plot_inv(diameters, age)
   @savefig inverse.png width=4in
   plt.show()

.. _ug-production-simulation:

Using the Production component in a Simulation
----------------------------------------------

In the above guides we have shown how to use the |Production| as a standalone object. However, the preferred way of initializing and using this object is as part of a |Simulation| object, where it is initialized in tandem with all other components in the Cratermaker project. This allows the components to be checked for self-consistency, and streamlines some of the initialization issues. Here we will demonstrate some typical use-cases for |Production|, including how the valid ranges are dynamically adjusted for the particular  |Surface| and how the Quasi-Monte Carlo functionality allows the user to mix real craters with known properties with randomly-generated craters into a simulation.

Emplace a population of craters drawn from a production function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The |sim.populate| is used to generate a population of craters drawn from whichever  |Production| component module is loaded into the simulation. You can pass either time interval over which to draw from the production function with either "age" or "time_start" and "time_end", or you can provide a pair of numbers in the form of (D,N) that represents the production function in the N(D) convention. For instance, to emplace a population of craters corresponding to the production function between 200 My to 100 My before present, you would do:

.. code-block:: python

    sim.populate(time_start=200, time_end=100)


The |sim.populate| method is really meant to be a helper method for the |sim.run| method, which is used to run a full simulation over a specified time interval. For instance, the |sim.populate| method does not track the time of emplacement or advance the time interval of the simulation.


Limiting the size range of the production population
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When used inside a |Simulation|, the upper and lower bounds on the sizes of craters returned by |production.function|


Quasi-Monte Carlo craters
-------------------------

For some applications, you may want to specify a predefined population of craters that should be emplaced at a certain time and location along with the random population of craters drawn from the production function. To accomplish this, you can initialize a |Simulation| object with an argument `quasimc_file`, which points to either a CSV or NetCDF file containing the information needed to emplace the craters. Alternatively, you  When this argumented is passed, the |sim.run| method will run in "Quasi-Monte Carlo mode". 

Quasi-Monte Carlo mode is a powerful tool that is very flexible. Here we will demonstrate the different ways you can specify craters in the input file and how they affect the behavior of the simulation with a simple donstration using several prominent lunar craters: South Pole-Aitken, Serenitatis, Nectaris, Crisium, Imbrium, Schrödinger, and Orientale. A version of this simulation type with 76 lunar basins and craters is included in the :ref:`gal-simulation`.

First, we we will create a CSV file called "qmc_inputs.csv" and populate it with columns indicating which arguments should be passed to the |Crater.maker| function. At a minimum, this requires at least one argument specifying size, such as :py:attr:`~cratermaker.components.crater.CraterFixed.diameter` or :py:attr:`~cratermaker.components.crater.CraterFixed.projectile_diameter` Usually you would include a location as well, which must be specified with "longitude" and "latitude" as separate columns, rather than the typical :py:attr:`~cratermaker.components.crater.CraterFixed.location` tuple that the method call would use. In addition, you typically would provide some indication of when in the simulation you want the crater to form. There are multiple ways to do this.


- |Crater.production_time|: This specifies the time when the crater should form along with an optional 1-σ standard deviation value. When passing this to |Crater.maker| this can be specified as either a single value or a pair of values, where the first value represents the time and the second the standard deviation, both in units of My. The actual time value will be drawn from a normal distribution. 

.. ipython:: python
    :okwarning:
    :suppress:

    from cratermaker import cleanup
    cleanup()


.. ipython:: python
   :okwarning:

   from cratermaker import Simulation

   sim = Simulation(gridlevel=5)
   qmc_list = []
   qmc_list.append(sim.Crater.maker(name="Imbrium", diameter=1321e3, location=(341.5, 37), production_time=(3922, 12)))
   sim.quasimc_craters = qmc_list # This will process the craters and give it a time values bsed on their production metadata
   print(sim.quasimc_craters[0])


You can indicate that you want a crater to form at a specific time, as defined by the chronology function of the Simulation's |Production| component. To do this, you would include the time in a column labeled "production_time". You may also indicate that you want to draw the crater's time from a distribution by supplying an aditional column labeled "production_time_stdev".  Here we will define dates in My before present that each of our lunar craters will form within the Neukum chronology function used for default simulation type. 

We will use the age of the Imbrium impact event reported by Nemchin et al. (2021) of 3922 My, however the true ages for the other craters are very poorly constrained. However, because the chronology function in |NPF| translates model ages into crater number densities, we can pick dates that place our craters in their correct relative sequence within the bombardment history of the early Moon, even if their true age is unknown. We begin our cratering simulation at 4310 My before present, and give that age to the largest and stratigraphically oldest lunar crater, the massive 2400 km diameter South Pole-Aitken basin. 

.. note::

    Always remember that crater chronology model ages are **not** true ages, especially for the period of time prior to 3900 My bp. It's best to think of the model age as a way of expressing crater number density (not the other way around), and is subject to change as better calibration data is obtained. Thus a model age of 4310 My bp just means N(20)=1000 craters per million sq. km. in the Neukum production function, and will likely be updated if and when samples of the South Pole-Aitken melt sheet are returned and dated.

Serenetitatis likely formed before Nectaris, which formed before Crisium and Imbrium, so we set ``production_time`` values for Serenitatis, Nectaris, and Crisium of 4220, 4170, and 4070 My, respectively. Both Schrödinger and Orientale post-date Imbrium, so we set their ``production_time`` values to 3860 and 3810 My, respectively. Therefore our 

.. ipython:: python
    :okwarning:
    :suppress:

    from pathlib import Path
    with Path.open("basins_exact_time.csv", "w") as f:
        f.write("name,latitude,longitude,diameter,production_time\n")
        f.write("South Pole-Aitken,-53,191,2400000,4310\n")
        f.write("Serenitatis,25.4,18.8,923000,4220\n")
        f.write("Nectaris,-15.6,35.1,885000,4170\n")
        f.write("Crisium,16.8,58.4,1076000,4070\n")
        f.write("Imbrium,37,341.5,1321000,3922\n")
        f.write("Schrödinger,-74.9,133.5,326000,3860\n")
        f.write("Orientale,-20.1,265.2,937000,3810\n")

.. ipython:: python
    :okwarning:

    from pathlib import Path

    with Path.open("basins_exact_time.csv", "r") as f:
        print(f.read())


Now we can pass this file into a Simulation, with a reduced resolution to make the example run more quickly. We will also cut off the crater population below 100 km by setting the |sim.smallest_crater| attribute, which will also speed up our simulation as the majority of craters will be small ones. We also pass a value to the rng_seed argument of the Simulation constructor so that the random population of craters drawn from the production function is the same each time this example is run.

.. code-block:: python

    from cratermaker import Simulation
    sim = Simulation(quasimc_file="basins_exact_time.csv",gridlevel=5, ask_overwrite=False, reset=True, rng_seed=298263286)
    sim.smallest_crater = 100e3
    sim.run(age=4310)
    sim.show3d(variable_name="face_elevation")

.. ipython:: python
   :okwarning:
   :suppress:

    from cratermaker import Simulation
    sim = Simulation(quasimc_file="basins_exact_time.csv",gridlevel=5, ask_overwrite=False, reset=True, rng_seed=298263286)
    sim.smallest_crater = 100e3
    sim.run(age=4310)


.. pyvista-plot::
   :caption: Quasi-Monte Carlo simulation with craters emplaced at specific times.
   :include-source: False 

    >>> # The pyvista-plot tool apparently doesn't recognize the show3d method, so we will just show how to do this manually with pyvista.
    >>> import pyvista as pv
    >>> import cratermaker as cm
    >>> sim = cm.Simulation(reset=False)
    >>> plotter = sim.pyvista_plotter(variable_name="face_elevation") 
    >>> plotter.show()



Emplacing a crater at a specific time range
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In some cases, we may be more confident that a crater forms within a specific time range. To do this, we can set the "production_time_stdev" and "production_time_high" columns in the input file, which will cause the Simulation to draw a random time for that crater in a normal distribution, where the mean is the midpoint between  


More Production examples
------------------------

See more examples at  :ref:`gal-production_and_montecarlo`

.. toctree::
   :maxdepth: 3
   :hidden:


References
----------

.. [#] Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55-86. https://doi.org/10.1023/A:1011989004263
.. [#] Ivanov, B.A., 2001. Mars/Moon Cratering Rate Ratio Estimates. *Space Science Reviews*, 96, 87-104. https://doi.org/10.1023/A:1011941121102
.. [#] Ivanov, B.A., Neukum, G., Wagner, R., 2001. Size-Frequency Distributions of Planetary Impact Craters and Asteroids, in: Collisional Processes in the Solar System. Springer Netherlands, Dordrecht, pp. 1-34. https://doi.org/10.1007/978-94-010-0712-2_1

.. ipython:: python
    :okwarning:
    :suppress:

    cleanup()