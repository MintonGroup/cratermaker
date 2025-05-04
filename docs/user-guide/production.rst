.. currentmodule:: cratermaker.production

.. _production:

##########
Production
##########

.. rubric:: by Austin Blevins and David A. Minton

Cratermaker's production module provides a robust way to compute the production function for craters and projectiles. The production function is defined as the cumulative number of craters greater than a given diameter per unit m\ :sup:`2` surface area.  There are currently two types of production funcitons that are available in the :ref:`api-Production` class: "neukum" and "powerlaw". The default production model depends on the value of the target. If the target is one of the inner planets (or Vesta or Ceres), the default will be "neukum". For all other bodies, it is "powerlaw".  The "neukum" production function comes in three different versions: "moon", "mars", or "projectile."  The "projectile" version is used unless the target is either the Moon or Mars. The default Cratermaker target body is the Moon, and therefore the default production function is "neukum" with the "moon" model. We can create a production function object and inspect its values using the following code:

.. ipython:: python
   :okwarning:

   from cratermaker import Production
   production = Production.maker()
   print(production)

The output shows that by default, the production function is the Neukum model for the Moon, which is described in Neukum, Ivanov, and Hartmann (2001) [1]_. It also reports that the "Generator type" is "crater." A production model can either a crater generator or a projectile generator, which means that the diameter values that are returned from the production function represent either the crater or projectile. The other crater-based model is the Mars version, which is based on the work of Ivanov (2001) [2]_. The projectile model is based on the work of Ivanov, Neukum, and Wagner (2001) [3]_.

The Neukum Moon model has a valid range of ages over which it is valid, which is also reported in the output. The valid range of ages is between 0 and 4.5 Ga, though caution must be used when interpreting ages older than about 3.9 Ga, as these are poorly calibrated. The total crater number density of the most ancient terrains on the Moon correspond to a model age of 4.31 Ga, though the actual crust could potentially be must older than that. Nevertheless, the Neukum model is a well-established and widely used model for inner solar system crater chronology.

The Production class has two primary methods that are used to compute the production function. The first is :func:`Production.function`, which computes the expected cumulative number density of craters greater than a given diameter of a surface of a given age. 
The second is :func:`Production.sample`, which samples crater diameters and ages from the production function using Monte Carlo methods.

Production function
===================
:func:`Production.function` returns the cumulative size-frequency distribution (CSFD) of craters over a given age range and crater diameter. It takes the following arguments:

- **diameter**: Crater diameter(s) in units of meters to compute corresponding cumulative number density value.
- **age**: Age in units of My relative to the present, used to compute the CSFD. Default is ``1.0``, corresponding to 1 Ma.
- **age_end**: ending age in units of My relative to the present, also used to compute the CSFD. Default is ``0.0``, which corresponds to the present day.

Example: Using Production.function
==================================

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
====================================================================================================


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

   @savefig ../_static/PowerLawExample.png width=4in
   ax.loglog(diameters, pf)


Example: Using Production.function_inverse
==========================================

:func:`Production.function_inverse` returns the age in My for a given crater number density and diameter. It has the following parameters:

- **diameter**: Diameter of the crater in m
- **cumulative_number_density**: Number density of craters per m\ :sup:`2` surface area greater than the input diameter.

For this example, we are going to use :func:`Production.function_inverse` to plot the age in Ma for a number density of 1e-6 craters per m² from 1 m to 1 km in diameter. 

.. ipython:: python
   :okwarning:

   from cratermaker import Production
   import matplotlib.pyplot as plt
   import numpy as np


   fig, ax = plt.subplots()
   ax.set_xlabel('D (m)')
   ax.set_ylabel('Age (Ma)')

   production = Production.maker("powerlaw",slope=-2.5)
   inv = production.function_inverse(diameter=np.arange(1,1000), cumulative_number_density=1e-6*np.ones(len(np.arange(1,1000))))

   plt.grid()
   plt.tight_layout()

   @savefig ../_static/inverse.png width=4in
   ax.loglog(np.arange(1,1000), inv)


Using the NeukumProduction class
================================

The :ref:`api-NeukumProduction` class can compute the Neukum production function for the Moon and Mars. The Neukum production function is  is computed using the model of Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86 for the Moon and Mars, with minor changes.

The main parameter is **model**, which is the specific model to use for the production function. Defaults to ``Moon``, but could also be ``Mars`` or ``Projectile``. The choice of ``model`` sets a number of parameters for the Neukum production. See :class:`cratermaker.NeukumProduction` for more, including references.

Example: Plot the Neukum projectile CSFD
========================================

.. ipython:: python
   :okwarning:

   from cratermaker import Production

   def plot_npf():
      import matplotlib.pyplot as plt
      import numpy as np
      import matplotlib.ticker as ticker

      fig, ax = plt.subplots(figsize=(4, 7))
      x_min = 1e-5
      x_max = 1e4
      y_min = 1e-7
      y_max = 1e13
      nD = 1000
      Dvals = np.logspace(np.log10(x_min), np.log10(x_max), num=nD)
      production = Production.maker("neukum", version="Projectile")
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
      inrange = (Dvals >= production.sfd_range[0]) & (Dvals <= production.sfd_range[1])
      lo = Dvals < production.sfd_range[0]
      hi = Dvals > production.sfd_range[1]
      t = 1.0
      Nvals = production.function(diameter=Dvals*1e3,age=t*1e3)
      Nvals *= 1e6 # convert from m^-2 to km^-2
      ax.plot(Dvals[inrange], Nvals[inrange], '-', color='black', linewidth=1.0, zorder=50)
      ax.plot(Dvals[lo], Nvals[lo], '-.', color='orange', linewidth=2.0, zorder=50)
      ax.plot(Dvals[hi], Nvals[hi], '-.', color='orange', linewidth=2.0, zorder=50)

      plt.tick_params(axis='y', which='minor')
      plt.tight_layout()

   @savefig ../_static/projectilecsfd.png width=4in
   plot_npf()



Example: Plot the NPF Chronology function for the Moon and Mars
===============================================================

.. ipython:: python
   :okwarning:

   from cratermaker import Production
   def plot_chronology_functions():
      import matplotlib.pyplot as plt
      import numpy as np

      fig, ax = plt.subplots(figsize=(8, 4))
      ax.set_yscale('log')
      ax.set_ylabel('$\\mathregular{N(1) (km^{-2})}$')
      ax.set_xlabel('Time (Gy ago)')
      ax.set_xlim(4.5, 0)
      moon = Production.maker("neukum", version="Moon")
      mars = Production.maker("neukum", version="Mars")
      tvals = np.linspace(4.5, 0.0, num=1000)
      N1_moon = moon.function(diameter=1000.0, age=tvals*1e3)*1e6
      N1_mars = mars.function(diameter=1000.0, age=tvals*1e3)*1e6
      ax.plot(tvals, N1_moon, '-', color='dimgrey', linewidth=2.0, zorder=50, label="Moon")
      ax.plot(tvals, N1_mars, '-', color='orange', linewidth=2.0, zorder=50, label="Mars")
      ax.legend()
      plt.tight_layout()
   
   @savefig ../_static/chronology.png width=5in
   plot_chronology_functions()

Example: Plot the NPF CSFD with isochrons at 1 Ma, 1 Ga, and 4 Ga
=================================================================

.. ipython:: python
   :okwarning:

   from cratermaker import Production

   def plot_isochrons():
      import matplotlib.pyplot as plt
      import numpy as np
      import matplotlib.ticker as ticker

      fig, axs = plt.subplots(1, 2, figsize=(8, 7))
      ax = {'Moon': axs[0], 'Mars': axs[1]}

      tvals = [0.01, 1.0, 4.0]
      x_min = 1e-3
      x_max = 1e4
      y_min = 1e-12
      y_max = 1e6
      nD = 1000
      Dvals = np.logspace(np.log10(x_min), np.log10(x_max), num=nD)
      for key in ax:
         production = Production.maker("neukum", version=key)
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
         inrange = (Dvals >= production.sfd_range[0]) & (Dvals <= production.sfd_range[1])
         lo = Dvals < production.sfd_range[0]
         hi = Dvals > production.sfd_range[1]
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

   @savefig ../_static/isochrons.png width=7in
   plot_isochrons()

:func:`Production.sample` allows you to sample crater diameters and ages from the production function (either power law or Neukum). This function can either sample from a given age range or from a given cumulative number/diameter pair, but not both. 

Example: Sample a power law and lunar Neukum Production Function
================================================================

.. ipython:: python
   :okwarning:

   from cratermaker import Production
   def sample_csfd():
      import matplotlib.pyplot as plt
      import numpy as np
      import matplotlib.ticker as ticker

      fig, axs = plt.subplots(1, 2, figsize=(8, 7))
      ax = {'Power Law': axs[0], 'NPF (Moon)': axs[1]}

      production = {
         'Power Law': Production.maker("powerlaw"),
         'NPF (Moon)': Production.maker("neukum", version="Moon")
      }

      from cratermaker import Target
      target = Target("Moon")
      area = 4 * np.pi * target.radius**2
      age = 4100.0
      x_min = 1e0
      x_max = 1e5
      y_min = 1e0
      y_max = 1e4
      diameter_range = (2e3, 10000e3) # Range of diameters to generate in m
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
            Nsampled = range(1, len(Dsampled)+1)
            ax[key].plot(Dsampled*1e-3, Nsampled, '-', color='cornflowerblue', linewidth=2.0, zorder=50, alpha=0.2)
         # Plot the production function
         Nvals = production[key].function(diameter=Dvals*1e3, age=age)
         Nvals *= area # convert from per unit area to total number
         ax[key].plot(Dvals, Nvals, '-', color='black', linewidth=3.0, zorder=50)

      plt.tick_params(axis='y', which='minor')
      plt.tight_layout()

   @savefig ../_static/samples.png width=7in
   sample_csfd()

.. toctree::
   :maxdepth: 2
   :hidden:

##########
References
##########


.. [1] Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. *Space Science Reviews*, 96, 55-86. https://doi.org/10.1023/A:1011989004263
.. [2] Ivanov, B.A., 2001. Mars/Moon Cratering Rate Ratio Estimates. *Space Science Reviews*, 96, 87-104. https://doi.org/10.1023/A:1011941121102
.. [3] Ivanov, B.A., Neukum, G., Wagner, R., 2001. Size-Frequency Distributions of Planetary Impact Craters and Asteroids. In *Collisional Processes in the Solar System*, Springer Netherlands, Dordrecht, pp. 1-34.  https://doi.org/10.1007/978-94-010-0712-2_1 