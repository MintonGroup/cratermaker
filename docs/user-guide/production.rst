.. currentmodule:: cratermaker.production

.. image:: ../_static/production_icon.svg
    :alt: Production
    :align: center
    :width: 600px


.. _ug-production:

Production
==========

Cratermaker's production module provides a robust way to compute the production function for craters and projectiles. The production function is defined as the cumulative number of craters greater than a given diameter per unit m² surface area.  There are currently two types of production funcitons that are available in the :ref:`Production <api-Production>` class: "neukum" and "powerlaw". The default production model depends on the value of the target. If the target is one of the inner planets (or Vesta or Ceres), the default will be "neukum". For all other bodies, it is "powerlaw".  The "neukum" production function comes in three different versions: "moon", "mars", or "projectile."  The "projectile" version is used unless the target is either the Moon or Mars. The default Cratermaker target body is the Moon, and therefore the default production function is "neukum" with the "moon" model. We can create a production function object and inspect its values using the following code:

.. ipython:: python
   :okwarning:

   from cratermaker import Production
   production = Production.maker()
   print(production)

The output shows that by default, the production function is the Neukum model for the Moon, which is described in Neukum, Ivanov, and Hartmann (2001) [#]_. It also reports that the "Generator type" is "crater." A production model can either a crater generator or a projectile generator, which means that the diameter values that are returned from the production function represent either the crater or projectile. The other crater-based model is the Mars version, which is based on the work of Ivanov (2001) [#]_. The projectile model is based on the work of Ivanov, Neukum, and Wagner (2001) [#]_.

The Neukum Moon model has a valid range of ages over which it is valid, which is also reported in the output. The valid range of ages is between 0 and 4.5 Ga, though caution must be used when interpreting ages older than about 3.9 Ga, as these are poorly calibrated. The total crater number density of the most ancient terrains on the Moon correspond to a model age of 4.31 Ga, though the actual crust could potentially be must older than that. Nevertheless, the Neukum model is a well-established and widely used model for inner solar system crater chronology.

The Production class has two primary methods that are used to compute the production function. The first is :func:`Production.function`, which computes the expected cumulative number density of craters greater than a given diameter of a surface of a given age. 
The second is :func:`Production.sample`, which samples crater diameters and ages from the production function using Monte Carlo methods.

Production function
-------------------

:func:`Production.function` returns the cumulative size-frequency distribution (CSFD) of craters over a given age range and crater diameter. It takes the following arguments:

- **diameter**: Crater diameter(s) in units of meters to compute corresponding cumulative number density value.
- **age**: Age in units of My relative to the present, used to compute the CSFD. Default is ``1.0``, corresponding to 1 Ma.
- **age_end**: ending age in units of My relative to the present, also used to compute the CSFD. Default is ``0.0``, which corresponds to the present day.

Example: Using Production.function
----------------------------------

For this example, we are going to use :func:`Production.function` to find the cumulative number of craters greater than 1 km² to form on a surface in 3 Gy using a power law with a slope of -2. 

.. ipython:: python
   :okwarning:

   from cratermaker import Production 
   production = Production.maker("powerlaw",slope=-2.0) 
   diameter = 1000.0 
   age = 3000.0  
   n = production.function(diameter=diameter,age=age)  
   print(f"{n*1e6:.2f} D>{diameter*1e-3:.0f} km craters per km² in {age*1e-3:.0f} Gy")

Example: Plot a power law production function with a slope of -3 for craters 1 m to 100 km over 1 Ga
----------------------------------------------------------------------------------------------------


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


Example: Using Production.function_inverse
------------------------------------------

:func:`Production.function_inverse` returns the age in My for a given crater number density and diameter. It has the following parameters:

- **diameter**: Diameter of the crater in m
- **cumulative_number_density**: Number density of craters per m² surface area greater than the input diameter.

For this example, we are going to use :func:`Production.function_inverse` to plot the age in Ma for a number density of 1e-6 craters per m² from 1 m to 1 km in diameter. 


.. ipython:: python
   :okwarning:

   from cratermaker import Production
   import matplotlib.pyplot as plt
   import numpy as np

   production = Production.maker("powerlaw", slope=-2.5)
   diameters = np.logspace(0, 3)
   cumulative_number_density = np.full_like(diameters, 1e-6)
   age = production.function_inverse(diameter=diameters, cumulative_number_density=cumulative_number_density)

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

More Production examples
------------------------

See more examples at  :ref:`gal-production`

.. toctree::
   :maxdepth: 2
   :hidden:


References
----------

.. [#] Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55-86. https://doi.org/10.1023/A:1011989004263
.. [#] Ivanov, B.A., 2001. Mars/Moon Cratering Rate Ratio Estimates. *Space Science Reviews*, 96, 87-104. https://doi.org/10.1023/A:1011941121102
.. [#] Ivanov, B.A., Neukum, G., Wagner, R., 2001. Size-Frequency Distributions of Planetary Impact Craters and Asteroids, in: Collisional Processes in the Solar System. Springer Netherlands, Dordrecht, pp. 1-34. https://doi.org/10.1007/978-94-010-0712-2_1