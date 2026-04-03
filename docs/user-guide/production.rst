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

In the above guides we have shown how to use the |Production| as a standalone object. However, the preferred way of initializing and using this object is as part of a |Simulation| object, where it is initialized in tandem with all other components in the Cratermaker project. This allows the components to be checked for self-consistency, and streamlines some of the initialization issues. Here we will demonstrate some typical use-cases for |Production|, including how the valid ranges are dynamically adjusted for the particular |Surface| and how the Quasi Monte Carlo functionality allows the user to mix real craters with known properties with randomly-generated craters into a simulation.

Emplace a population of craters drawn from a production function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The |sim.populate| is used to generate a population of craters drawn from whichever |Production| component module is loaded into the simulation. You can pass either time interval over which to draw from the production function with either "age" or "time_start" and "time_end", or you can provide a pair of numbers in the form of (D,N) that represents the production function in the N(D) convention. For instance, to emplace a population of craters corresponding to the production function between 200 My to 100 My before present, you would do:

.. code-block:: python

    sim.populate(time_start=200, time_end=100)


The |sim.populate| method is really meant to be a helper method for the |sim.run| method, which is used to run a full simulation over a specified time interval. For instance, the |sim.populate| method does not track the time of emplacement or advance the time interval of the simulation.


Limiting the size range of the production population
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When used inside a |Simulation|, the upper and lower bounds on the sizes of craters returned by |production.function| are set by the resolution of the associated |Surface| object. The smallest crater that is directly modeled (outside of any subpixel models) is the size of the local face. The largest size is the diameter of the target. These can be adjusted by the user and are access by the |sim.smallest_crater| and |sim.largest_crater| attributes.


.. ipython:: python
    :okwarning:
    :suppress:

    from cratermaker import cleanup
    cleanup()

.. ipython:: python
   :okwarning:


   from cratermaker import Simulation

   sim = Simulation(gridlevel=5)
   sim.smallest_crater
   sim.largest_crater

.. _ug-production-quasimc:

Quasi Monte Carlo craters
-------------------------

For some applications, you may want to specify a predefined population of craters that should be emplaced at a certain time and location along with the random population of craters drawn from the production function. To accomplish this, you can initialize a |Simulation| object with an argument `quasimc_file`, which points to either a CSV or NetCDF file containing the information needed to emplace the craters. Alternatively, you  When this argumented is passed, the |sim.run| method will run in "Quasi Monte Carlo mode". 

Quasi Monte Carlo mode is a powerful tool that is very flexible. Here we will demonstrate the different ways you can specify craters in either a file via the |production.quasimc_file| property, or a list of |Crater| objects provided to the |production.quasimc_craters| property. The |Simulation| class also has its own wrappers for these two properties, |sim.quasimc_file| and |sim.quasimc_craters| for convenience.  

Here we have some demonstrations of the different ways to set Quasi Monte Carlo craters, and how different parameters affect the behavior of the simulation with a simple donstration using several prominent lunar craters: South Pole-Aitken, Serenitatis, Nectaris, Crisium, Imbrium, Schrödinger, and Orientale. A version of this simulation type with 76 lunar basins and craters is included in the :ref:`gal-simulation`.

First, we we will create a CSV file called "qmc_inputs.csv" and populate it with columns indicating which arguments should be passed to the |Crater.maker| function. At a minimum, this requires at least one argument specifying size, such as :py:attr:`~cratermaker.components.crater.CraterFixed.diameter` or :py:attr:`~cratermaker.components.crater.CraterFixed.projectile_diameter` Usually you would include a location as well, which must be specified with "longitude" and "latitude" as separate columns, rather than the typical :py:attr:`~cratermaker.components.crater.CraterFixed.location` tuple that the method call would use. In addition, you typically would provide some indication of when in the simulation you want the crater to form. There are multiple ways to do this, which we will demonstrate using ancient lunar basins.

.. note::

   We used the following sources for constructing our Quasi Monte Carlo constraints. This is an updated set of constraints from those used by Blevins et al. (2025) [#]_.

   - The oldest crater in the file is South Pole Aitken basin, which we have given an N(20) of 999.8 per 10⁶ km² and a sequence number of 0. This puts it just after our chosen simulation start time in terms of NPF model years.
   - The 74 largest and oldest crater names, sizes, and impact locations are those determined using GRAIL-derived gravity data from Neumann et al. (2015) [#]_
   - Sequence numbers are determined by stratigraphic relationships from Wilhelms et al. (1987) [#]_ and Fassett et al. (2012) [#]_
   - Fecunditatis and Australe North are assigned N(90) values from Evans et al. (2018) [#]_.
   - We do not consider the evidence for Procellarum being an impact basin as particularly strong, and so we do not include it.
   - Basins with N(20) values or those with stratigraphic sequences but no ages are taken from Orgel et al. (2018) [#]_.
   - Imbrium is set to 3922±12 My based on U-Pb ages of Apollo 14 and 15 impact breccias from Nemchin et al. (2021) [#]_.
   - The age of Copernicus of 800±15 My is from Bogard et al. (1994) [#]_, and is based on :sup:`39`\ Ar - :sup:`40`\ Ar ages of Apollo 12 soil samples that are thought to have been reset by ejecta from Copernicus.
   - The age of Tycho of 109±4 My is based on cosmic ray exposure dating of the South Massif landslide deposit in Taurus Littrow Valley from Apollo 17 by Drozd et al. (1977) [#]_.
   - Serenetatis is given an age of 4250 My based on modeling the origin of the Apollo 17 troctolite 76535 (which we think is The Most Interesting Rock from the Moon) by Bjonnes et al. (2025) [#]_. However, this is still a proposed age and because it would represent a date that is not part of the calibration used to develop the Neukum chronology function, then it may not be consistent with its position in the lunar chronology. Alternatively, Orgel et al.'s N(20)=334±73 could be used instead.
   


Specifying an emplacement time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Setting |crater.production_time| specifies the time when the crater should form along with an optional 1σ standard deviation value. When passing this to |Crater.maker| this can be specified as either a single value or a pair of values, where the first value represents the time and the second the standard deviation, both in units of My. The actual time value will be drawn from a normal distribution, unless you only pass a single value to |crater.production_time|.  

.. ipython:: python
    :okwarning:
    :suppress:

    from cratermaker import cleanup
    cleanup()


.. ipython:: python
   :okwarning:

   from cratermaker import Simulation

   sim = Simulation(gridlevel=5)
   sim.quasimc_craters.append(sim.Crater.maker(name="Imbrium", diameter=1321e3, location=(341.5, 37), production_time=(3922, 12)))
   print(sim.quasimc_craters[-1])


The equivalent in a CSV file would be:

.. csv-table:: qmc_imbrium.csv
   :header: "name","latitude","longitude","diameter","production_time","production_time_stdev"
   :widths: auto

   Imbrium,37,341.5,1321000,3922,12


.. note::

    Always remember that crater chronology model ages are **not** true ages, especially for the period of time prior to 3900 My bp. It's best to think of the model age as a way of expressing crater number density (not the other way around), and is subject to change as better calibration data is obtained. Thus a model age of 4310 My bp just means N(20)=999.8 craters per million sq. km. in the Neukum chronology function, which in reality could represent anything from the high 4400s to the mid 4200s. Getting a more precise date for the formation of South Pole-Aitken basin, which is stratigraphically the oldest surface feature known on the Moon, will likely have to wait until samples of its melt sheet are returned and dated.

Specifying a production crater number density
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rather than a specific time, you can specify a crater's "age" by its position on the production function using N(D) convention. This is done by setting the |crater.production_ND| attribute, which is a tuple of (D, N, N_stdev), where D is the reference diameter in km, N is the cumulative number density of craters per 10⁶ km² greater than D, and N_stdev is the 1σ standard deviation of that number density. Here we will use the N(20) value of the Nectaris basin obtained by Orgel et al. (2020) using the Buffered Non-Sparseness Correction technique, which yielded an N(20)=172±20 per 10⁶ km².

.. note::

   The N(D) convention takes D in units of km, rather than the units of m used when defining craters. To see what unit system is current set, see the |production.N_D_units|. Other conventions can be set changing |production.D_conversion_factor| and |production.N_conversion_factor|.


.. ipython:: python
   :okwarning:


   sim.quasimc_craters.append(sim.Crater.maker(name="Nectaris", diameter=885e3, location=(35.1, -15.6), production_ND=(20, 172, 20)))
   print(sim.quasimc_craters[0]) # The quasimc_craters list is sorted by age, from oldest to youngest. So even though Nectaris was added after Imbrium, it is actually the oldest crater in the list, so it is at index 0.


The equivalent in a CSV file would be:

.. csv-table:: qmc_with_nectaris.csv
   :header: "name","latitude","longitude","diameter","production_time","production_time_stdev","production_D","production_N","production_N_stdev"
   :widths: auto

   Nectaris,-15.6,35.1,885000,,,20,172,20
   Imbrium,37,341.5,1321000,3922,12,,,

Specifying a sequence constraint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For many craters, it is difficult to constrain their formation time, however we can bracket them stratigraphically relative to other craters. To accomplish this in Cratermaker, we can supply an integer value to the |crater.production_sequence| attribute. Craters that have a sequence number must form after craters with a lower sequence number. A sequence number can be added to any crater, even if it has another piece of production metadata such as |crater.production_time| or |crater.production_ND|, or has none. The only constraint is that at least one crater with the lowest sequence number must have either |crater.production_time| or |crater.production_ND| set, otherwise there is no way to determine the earliest time boundary of the sequences.

In this example, will add the basins Fecunditatis and Vaporum. Fecunditatis was determined by Evans et al. (2018) to have  N(90)=10±4 per 10⁶ km². Vaporum is stratigraphically older than Fecunditatis, but younger than South Pole-Aitken (SPA). SPA is stratigraphically the oldest crater on the Moon, so we give it a value of |crater.production_sequence| of 0. Because it is the oldest, we need an additional constraint. Its absolute formation age is currently not well known. It must be younger than the Moon, so is likely less than 4500 My. There are materials older than about 4250 My that are attributed to basins other than SPA, so it is likely older than that. We will therefore constrain it to an N(20)=999.8 per 10⁶ km², which is the value of the production function at just after 4310 My bp, which we use as the model starting time of our simulation. Fecunditatis is younger than Vaporum, so we give it a sequence number of 10, and Vaporum gets a sequence number of 5, and no other constraints. We can also add sequence constraints to Nectaris and Imbrium.

Our input file for this set of craters would look like this:

.. ipython:: python
    :okwarning:
    :suppress:

    from pathlib import Path
    with Path.open("qmc_selected_basins.csv", "w") as f:
        f.write("name,latitude,longitude,diameter,production_time,production_time_stdev,production_D,production_N,production_N_stdev,production_sequence\n")
        f.write("South Pole-Aitken,-53,191,2400000,,,20,999,,0\n")
        f.write("Vaporum,14.2,3.1,410000,,,,,,5\n")
        f.write("Fecunditatis,-4.6,52,690000,,,90,10,4,10\n")
        f.write("Nectaris,-15.6,35.1,885000,,,20,172,20,40\n")
        f.write("Imbrium,37,341.5,1321000,3922,12,,,,100\n")


.. csv-table:: qmc_selected_basins.csv
   :file: ../qmc_selected_basins.csv
   :widths: auto


.. ipython:: python
   :okwarning:

   sim = Simulation(gridlevel=5, quasimc_file="qmc_selected_basins.csv", ask_overwrite=False)
   for crater in sim.quasimc_craters:
       N20 = sim.production.function(diameter=20e3, age=crater.time) * 1e12
       N64 = sim.production.function(diameter=64e3, age=crater.time) * 1e12
       n20text = f"N(20) = {N20:.1f}"
       n64text = f"N(64) = {N64:.1f}"
       stext = f"{crater.production_sequence}" if crater.production_sequence is not None else "-"
       print(f"{crater.name:22}: {crater.time:.1f} My {n20text:15} {n64text:14} {stext:5}")

.. note::

   The numerical value of the |crater.production_sequence| is only important if there are no other production metadata constraints, such as in the Vaporum example. In that case, the N(1) range corresponding to the sequence is an interpolation between bordering sequences. So in the above case, Vaporum's nominal "age" (in N(1) value) would be determined by extrapolating the N(1) values of craters with sequence 10 (Fecunditatis) and sequence 0 (SPA). With a value of 5, Vaporum would have a nominal N(1) value exactly in the middle, though its actual value will still be drawn from a normal distribution.


Adjusting the automatic crater size limit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, when a list of craters with production metadata is loaded into |production.quasimc_craters|, the upper range of the randomly-generated craters is adjusted so that it is limited to the size of the smallest crater in the |production.quasimc_craters| list. However, sometimes this is not desirable. For instance, suppose we want to add Copernicus to our list. If we didn't change the range of craters, this would mean that the simulation would never create random craters larger than the 96 km diameter of Copernicus. We can therefore adjust the |sim.largest_crater| attribute as needed. 




.. ipython:: python
    :okwarning:
    :suppress:

    from pathlib import Path
    with Path.open("qmc_with_copernicus.csv", "w") as f:
        f.write("name,latitude,longitude,diameter,production_time,production_time_stdev,production_D,production_N,production_N_stdev,production_sequence\n")
        f.write("South Pole-Aitken,-53,191,2400000,,,20,999,,0\n")
        f.write("Vaporum,14.2,3.1,410000,,,,,,5\n")
        f.write("Fecunditatis,-4.6,52,690000,,,90,10,4,10\n")
        f.write("Nectaris,-15.6,35.1,885000,,,20,172,20,40\n")
        f.write("Imbrium,37,341.5,1321000,3922,12,,,,100\n")   
        f.write("Copernicus,9.6209,339.9214,96070,800.0,15,,,,200\n")   


.. csv-table:: qmc_with_copernicus.csv
   :file: ../qmc_with_copernicus.csv
   :widths: auto


.. ipython:: python
   :okwarning:

   sim.quasimc_file="qmc_with_copernicus.csv"
   print(f"Largest random crater: {sim.largest_crater * 1e-3:.1f} km")

   # Now set the largest crater 
   sim.largest_crater = 200e3

   print(f"Largest random crater: {sim.largest_crater * 1e-3:.1f} km")




Modifying the Quasi Monte Carlo crater list
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The list of craters stored in |production.quasimc_craters| can be modified at any time. Reprocessing is triggered any time a crater in the list is found without a |crater.time| value. For instance, if you append to the list, the new crater will be added and the list will be reprocessed to extract new |crater.time| values.  Note that all craters in the list with production metadata will be reprocessed to generate a new |crater.time| value, even if they had one previously. 

.. ipython:: python
   :okwarning:

   # Append the crater Tycho to the list

   # Before adding Tycho
   for c in sim.quasimc_craters:
       print(f"{c.name}: Emplacement time {c.time} My bp")

   # Append Tycho to the list
   sim.quasimc_craters.append(
       sim.Crater.maker(
           name="Tycho",
           location=(348.7847, -43.2958),
           diameter=85294,
           production_time=(109, 4.0),
           production_sequence=300,
       )
   )

   # After adding Tycho
   for c in sim.quasimc_craters:
       print(f"{c.name}: Emplacement time {c.time} My bp")


In addition, setting a new file to |production.quasimc_file| will replace the existing |production.quasimc_craters| list and reprocess the new one.




Merging the Quasi Monte Carlo craters with random craters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When running in Quasi Monte Carlo mode, the randomly-generated and Quasi Monte Carlo generated craters are combined using the |production.quasimc_merge| (or its convenience wrapper |sim.quasimc_merge|). The arguments to this function are "craters", which holds a list of |Crater| objects that you have previously created, and either "time_start" and "time_end" or "N_D" and "N_D_end". "N_D" and "N_D_end" are the N(D) values expressed as a tuple of (D, N), where D is in km and N is in units of per 10⁶ km² (unless these have been modified by setting |production.D_conversion_factor| and/or |production.N_conversion_factor|).  Internally, the "N_D" and "N_D_end" values are converted to "time_start" and "time_end" values using the |production.age_from_D_N| function.

The |production.quasimc_merge| method will do a number of things to process both the input list of |Crater| objects given by the "craters" argument as well as the stored |production.quasimc_craters| list to produce a consistent merge. First, it checks all of the "craters" list to see if they have a |crater.time| value set, set them using the |production.compute_time| method if they don't have one or drop them if the pre-existing |crater.time| value falls outside the range. Then it successively checks the two lists for overlapping diameters, starting with the smallest craters in both. If a crater from the input (i.e. random) list has a diameter greater than the smallest of the quasimc list, then the quasimc crater takes its place, and then the next largest quasimc crater is checked against the next largest input crater. This continues intil both lists are exhausted.  The purpose of this is to ensure that the total number of craters produced over a time interval with a |production.quasimc_craters| component is roughly consistent with the total number that would be produced without. 

To demonstrate this functionality, consider the case where we want to model the total number of craters produced during the Copernican period, but we also want to include Copernicus and Tycho. First let's set up a low resolution |Simulation| object and not specify any Quasi Monte Carlo craters, but give it a random seed for repeatability

.. ipython:: python
    :okwarning:
    :suppress:

    from cratermaker import cleanup
    cleanup()


.. ipython:: python
   :okwarning:

   from cratermaker import Simulation

   sim = Simulation(gridlevel=5, rng_seed=11123537)


Next we will use |sim.populate| to generate a set of craters for the last 1 billion years. This method will return the list of randomly generated emplaced craters, which we will sort by diameter and plot them

.. ipython:: python
   :okwarning:

   random_craters = sim.populate(time_start=1000)
   random_craters.sort(key=lambda c: c.diameter, reverse=True)
   print()
   print(f"Random crater count: {len(random_craters)}")
   name = "Random"
   for c in random_craters:
       print(f"{name:10} crater: diameter = {c.diameter * 1e-3:.2f} km, time = {c.time:.2f} Myr")


Now we will add Copernicus and Tycho to the |sim.quasimc_craters| list.

.. ipython:: python
   :okwarning:
   
   sim.quasimc_craters = [
       sim.Crater.maker(name="Copernicus", location=(339.9214, 9.6209), diameter=96070, production_time=(800.0, 15), production_sequence=200),
       sim.Crater.maker(name="Tycho", location=(348.7847, -43.2958), diameter=85294, production_time=(109, 4.0), production_sequence=300)
   ]


We then merge the two lists and print the resulting list.

.. ipython:: python
   :okwarning:

   merged_craters = sim.quasimc_merge(craters=random_craters, time_start=1000)
   
   print(f"Merged crater count: {len(merged_craters)}")
   merged_craters.sort(key=lambda c: c.diameter, reverse=True)
   for c in merged_craters:
       if c.name is not None:
           name = c.name
       else:
           name = "Random"
       print(f"{name:10}: diameter = {c.diameter * 1e-3:.2f} km, time = {c.time:.2f} Myr")


Notice that the total number of craters in the merged list is the same as the original random list, but that some of the random craters have been substituted for quasimc ones. 

.. note::
   Usually you would not need to call the |production.quasimc_merge| or |sim.quasimc_merge| method directly. This method is automatically called by |sim.run| and |sim.populate| when creating the population of craters to emplace during a running simulation. 

   

More Production examples
------------------------

See more examples at  :ref:`gal-production_and_montecarlo` and :ref:`gal-simulation`.

.. toctree::
   :maxdepth: 3
   :hidden:


References
----------

.. [#] Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55-86. `doi:10.1023/A:1011989004263 <https://doi.org/10.1023/A:1011989004263>`_
.. [#] Ivanov, B.A., 2001. Mars/Moon Cratering Rate Ratio Estimates. *Space Science Reviews*, 96, 87-104. `doi:10.1023/A:1011941121102 <https://doi.org/10.1023/A:1011941121102>`_
.. [#] Ivanov, B.A., Neukum, G., Wagner, R., 2001. Size-Frequency Distributions of Planetary Impact Craters and Asteroids, in: Collisional Processes in the Solar System. Springer Netherlands, Dordrecht, pp. 1-34. `doi:10.1007/978-94-010-0712-2_1 <https://doi.org/10.1007/978-94-010-0712-2_1>`_
.. [#] Blevins, A.M., Minton, D.A., Huang, Y.H., Du, J., Tremblay, M.M., Fassett, C.I., 2025. Constraining the Source Craters of Apollo Impact Melts. JGR Planets 130, e2025JE009137. `doi:10.1029/2025JE009137 <https://doi.org/10.1029/2025JE009137>`_
.. [#] Neumann, G.A., Zuber, M.T., Wieczorek, M.A., Head, J.W., Baker, D.M.H., Solomon, S.C., Smith, D.E., Lemoine, F.G., Mazarico, E., Sabaka, T.J., Goossens, S.J., Melosh, H.J., Phillips, R.J., Asmar, S.W., Konopliv, A.S., Williams, J.G., Sori, M.M., Soderblom, J.M., Miljković, K., Andrews-Hanna, J.C., Nimmo, F., Kiefer, W.S., 2015. Lunar impact basins revealed by Gravity Recovery and Interior Laboratory measurements. Science Advances 1, e1500852–e1500852. `doi: 10.1126/sciadv.1500852 <https://doi.org/10.1126/sciadv.1500852>`_
.. [#] Wilhelms, D.E., McCauley, J.F., Trask, N.J., 1987. The geologic history of the moon. Washington : U.S. G.P.O. ; Denver.
.. [#] Fassett, C.I., Head, J.W., Kadish, S.J., Mazarico, E., Neumann, G.A., Smith, D.E., Zuber, M.T., 2012. Lunar impact basins: Stratigraphy, sequence and ages from superposed impact crater populations measured from Lunar Orbiter Laser Altimeter (LOLA) data. JGR 117, 1–13. `doi:10.1029/2011JE003951 <https://doi.org/10.1029/2011JE003951>`_
.. [#] Evans, A.J., Hanna, J.C.A., Head, J.W., Soderblom, J.M., Solomon, S.C., Zuber, M.T., 2018. Reexamination of Early Lunar Chronology With GRAIL Data: Terranes, Basins, and Impact Fluxes. J. Geophys. Res. 5, 313–1617. `doi:10.1029/2017JE005421 <https://doi.org/10.1029/2017JE005421>`_
.. [#] Orgel, C., Michael, G., Fassett, C.I., van der Bogert, C.H., Riedel, C., Kneissl, T., Hiesinger, H., 2018. Ancient Bombardment of the Inner Solar System: Reinvestigation of the “Fingerprints” of Different Impactor Populations on the Lunar Surface. J. Geophys. Res. 123, 748–762. `doi: 10.1002/2017JE005451 <https://doi.org/10.1002/2017JE005451>`_
.. [#] Nemchin, A.A., Long, T., Jolliff, B.L., Wan, Y., Snape, J.F., Zeigler, R., Grange, M.L., Liu, D., Whitehouse, M.J., Timms, N.E., Jourdan, F., 2021. Ages of lunar impact breccias: Limits for timing of the Imbrium impact. Geochemistry 81, 125683. `doi:10.1016/j.chemer.2020.125683 <https://doi.org/10.1016/j.chemer.2020.125683>`_
.. [#] Bogard, D.D., Garrison, D.H., Shih, C.Y., Nyquist, L.E., 1994. 39Ar-40Ar dating of two lunar granites: The age of Copernicus. Geochimica et Cosmochimica Acta 58, 3093–3100. `doi:10.1016/0016-7037(94)90181-3 <https://doi.org/10.1016/0016-7037(94)90181-3>`_
.. [#] Drozd, R.J., Hohenberg, C.M., Morgan, C.J., Podosek, F.A., Wroge, M.L., 1977. Cosmic-ray exposure history at Taurus-Littrow, in: Proceedings of the Lunar Science Conference. pp. 3027–3043. `abstract <https://ui.adsabs.harvard.edu/abs/1977LPSC....8.3027D/abstract>`_
.. [#] Bjonnes, E., Johnson, B.C., Broquet, A., Garrick‐Bethell, I., Andrews‐Hanna, J.C., Wakita, S., Kiefer, W.S., 2025. Evidence for an Early Formation of Serenitatis Basin at 4.25 Ga Shifts Lunar Chronology. Geophysical Research Letters 52, e2025GL116654. `doi:10.1029/2025GL116654 <https://doi.org/10.1029/2025GL116654>`_



.. ipython:: python
    :okwarning:
    :suppress:

    cleanup()