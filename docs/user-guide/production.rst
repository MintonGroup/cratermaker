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

Quasi-Monte Carlo craters
-------------------------

For some applications, you may want to specify a predefined population of craters that should be emplaced at a certain time and location along with the random population of craters drawn from the production function. To accomplish this, you can initialize a |Simulation| object with an argument `quasimc_file`, which points to either a CSV or NetCDF file containing the information needed to emplace the craters. Alternatively, you  When this argumented is passed, the |sim.run| method will run in "Quasi-Monte Carlo mode". 

Quasi-Monte Carlo mode is a powerful tool that is very flexible. Here we will demonstrate the different ways you can specify craters in the input file and how they affect the behavior of the simulation with a simple donstration using several prominent lunar craters: South Pole-Aitken, Serenitatis, Nectaris, Crisium, Imbrium, Schrödinger, and Orientale. A version of this simulation type with 76 lunar basins and craters is included in the :ref:`gal-simulation`.

First, we we will create a CSV file called "qmc_inputs.csv" and populate it with columns indicating which arguments should be passed to the |Crater.maker| function. At a minimum, this requires at least one argument specifying size, such as :py:attr:`~cratermaker.components.crater.CraterFixed.diameter` or :py:attr:`~cratermaker.components.crater.CraterFixed.projectile_diameter` Usually you would include a location as well, which must be specified with "longitude" and "latitude" as separate columns, rather than the typical :py:attr:`~cratermaker.components.crater.CraterFixed.location` tuple that the method call would use. In addition, you typically would provide some indication of when in the simulation you want the crater to form. There are multiple ways to do this.


Specifying an emplacement time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Setting |Crater.production_time| specifies the time when the crater should form along with an optional 1σ standard deviation value. When passing this to |Crater.maker| this can be specified as either a single value or a pair of values, where the first value represents the time and the second the standard deviation, both in units of My. The actual time value will be drawn from a normal distribution, unless you only pass a single value to |Crater.production_time|.  We will use the age of the Imbrium impact event reported by Nemchin et al. (2021) [#]_ of 3922±12 My.

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
   print(sim.quasimc_craters[-1])


The equivalent in a CSV file would be:

.. csv-table:: qmc_imbrium
   :header: "name","latitude","longitude","diameter","production_time","production_time_stdev"
   :widths: auto

   Imbrium,37,341.5,1321000,3922,12


.. note::

    Always remember that crater chronology model ages are **not** true ages, especially for the period of time prior to 3900 My bp. It's best to think of the model age as a way of expressing crater number density (not the other way around), and is subject to change as better calibration data is obtained. Thus a model age of 4310 My bp just means N(20)=999.8 craters per million sq. km. in the Neukum chronology function, which in reality could represent anything from the high 4400s to the mid 4200s. The formation age of the ancient South Pole-Aitken basins, which is stratigraphically the oldest surface feature known on the Moon, has to wait until samples of its melt sheet are returned and dated.

Specifying a production crater number density
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rather than a specific time, you can specify a crater's "age" by its position on the production function using N(D) convention. This is done by setting the |Crater.production_ND| attribute, which is a tuple of (D, N, N_stdev), where D is the reference diameter in km, N is the cumulative number density of craters per 10⁶ km² greater than D, and N_stdev is the 1σ standard deviation of that number density. Here we will use the N(20) value of the Nectaris basin obtained by Orgel et al. (2020) [#]_ using the Buffered Non-Sparseness Correction technique, which yielded an N(20)=172±20 per 10⁶ km².

.. note::

   The N(D) convention takes D in units of km, rather than the units of m used when defining craters. To see what unit system is current set, see the |production.N_D_units|. Other conventions can be set changing |production.D_conversion_factor| and |production.N_conversion_factor|.


.. ipython:: python
   :okwarning:


   qmc_list.append(sim.Crater.maker(name="Nectaris", diameter=885e3, location=(35.1, -15.6), production_ND=(20, 172, 20)))
   sim.quasimc_craters = qmc_list # This will process the craters and give it a time values bsed on their production metadata
   print(sim.quasimc_craters[-1])


The equivalent in a CSV file would be:

.. csv-table:: qmc_nectaris
   :header: "name","latitude","longitude","diameter","production_time","production_time_stdev","production_D","production_N","production_N_stdev"
   :widths: auto

   Nectaris,-15.6,35.1,885000,,,20,172,20
   Imbrium,37,341.5,1321000,3922,12,,,

Specifying a sequence constraint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For many craters, it is difficult to constrain their formation time, however we can bracket them stratigraphically relative to other craters. To accomplish this in Cratermaker, we can supply an integer value to the |Crater.production_sequence| attribute. Craters that have a sequence number must form after craters with a lower sequence number. A sequence number can be added to any crater, even if it has another piece of production metadata such as |Crater.production_time| or |Crater.production_ND|, or has none. The only constraint is that at least one crater with the lowest sequence number must have either |Crater.production_time| or |Crater.production_ND| set, otherwise there is no way to determine the earliest time boundary of the sequences.

In this example, will add the basins Fecunditatis and Vaporum. Fecunditatis was determined by Orgel et al. (2020) to have  N(90)=10±4 per 10⁶ km². Vaporum is stratigraphically older than Fecunditatis, but younger than South Pole-Aitken (SPA). SPA is stratigraphically the oldest crater on the Moon, so we give it a value of |Crater.production_sequence| of 0. Because it is the oldest, we need an additional constraint. Its absolute formation age is currently not well known. It must be younger than the Moon, so is likely less than 4500 My. There are materials older than about 4250 My that are attributed to basins other than SPA, so it is likely older than that. We will therefore constrain it to an N(20)=999 per 10⁶ km², which is the value of the production function at 4310 My bp. Fecunditatis is younger than Vaporum, so we give it a sequence number of 10, and Vaporum gets a sequence number of 5, and no other constraints. We can also add sequence constraints to Nectaris and Imbrium.

Our input file for this set of craters would look like this:

.. csv-table:: qmc_selected_basins
   :header: "name","latitude","longitude","diameter","production_time","production_time_stdev","production_D","production_N","production_N_stdev","production_sequence"
   :widths: auto

   South Pole-Aitken,-53,191,2400000,,,20,999,,0
   Vaporum,14.2,3.1,410000,,,,,,5
   Fecunditatis,-4.6,52,690000,,,90,10,4,10
   Nectaris,-15.6,35.1,885000,,,20,172,20,40
   Imbrium,37,341.5,1321000,3922,12,,,,100



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

   The numerical value of the |Crater.production_sequence| is only important if there are no other production metadata constraints, such as in the Vaporum example. In that case, the N(1) range corresponding to the sequence is an interpolation between bordering sequences. So in the above case, Vaporum's nominal "age" (in N(1) value) would be determined by extrapolating the N(1) values of craters with sequence 10 (Fecunditatis) and sequence 0 (SPA). With a value of 5, Vaporum would have a nominal N(1) value exactly in the middle, though its actual value will still be drawn from a normal distribution.


Adjusting the automatic crater size limit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, when a list of craters with production metadata is loaded into |production.quasimc_craters|, the upper range of the randomly-generated craters is adjusted so that it is limited to the size of the smallest crater in the |production.quasimc_craters| list. However, sometimes this is not desirable. For instance, suppose we want to add Copernicus to our list. If we didn't change the range of craters, this would mean that the simulation would never create random craters larger than the 96 km diameter of Copernicus. We can therefore adjust the |sim.largest_crater| attribute as needed. 


.. csv-table:: qmc_with_copernicus
   :header: "name","latitude","longitude","diameter","production_time","production_time_stdev","production_D","production_N","production_N_stdev","production_sequence"
   :widths: auto

   South Pole-Aitken,-53,191,2400000,,,20,999,,0
   Vaporum,14.2,3.1,410000,,,,,,5
   Fecunditatis,-4.6,52,690000,,,90,10,4,10
   Nectaris,-15.6,35.1,885000,,,20,172,20,40
   Imbrium,37,341.5,1321000,3922,12,,,,100
   Copernicus,9.6209,339.9214,96070,800.0,15,,,,200


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


.. ipython:: python
   :okwarning:

   sim.quasimc_file="qmc_with_copernicus.csv"
   print(f"Largest random crater: {sim.largest_crater * 1e-3:.1f} km")

   # Now set the largest crater 
   sim.largest_crater = 200e3

   print(f"Largest random crater: {sim.largest_crater * 1e-3:.1f} km")


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
.. [#] Nemchin, A.A., Long, T., Jolliff, B.L., Wan, Y., Snape, J.F., Zeigler, R., Grange, M.L., Liu, D., Whitehouse, M.J., Timms, N.E., Jourdan, F., 2021. Ages of lunar impact breccias: Limits for timing of the Imbrium impact. Geochemistry 81, 125683. `doi:10.1016/j.chemer.2020.125683 <https://doi.org/10.1016/j.chemer.2020.125683>`_
.. [#] Orgel, C., Michael, G., Fassett, C.I., van der Bogert, C.H., Riedel, C., Kneissl, T., Hiesinger, H., 2018. Ancient Bombardment of the Inner Solar System: Reinvestigation of the “Fingerprints” of Different Impactor Populations on the Lunar Surface. J. Geophys. Res. 123, 748–762. `doi:10.1002/2017JE005451 <https://doi.org/10.1002/2017JE005451>`_



.. ipython:: python
    :okwarning:
    :suppress:

    cleanup()