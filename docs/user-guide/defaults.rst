Default Behavior
================
As Cratermaker is designed to be easy to use, all of its component classes are built to be invoked without any arguments. When a ``Simulation`` object is created as in the examples above, a set of component objects are created with their default parameters. It is important to understand what the default behavior of these classes is, as this will affect the results of your simulation.

- **Target**: The default target is the Moon. There are a number of known bodies that can be selected as targets, including some jovian satellites, and small bodies, but currently few of them have been tested. 
- **Scaling**: The default scaling model is called ``default``, as it is similar to the one used by the Cratered Terrain Evolution Model (CTEM) that is a progenitor to Cratermaker [richardson2009]_. The projectile to transient scaling model is mostly based on Holsapple (1993) [holsapple1993]_ with some additional scaling parameters for ice given by Kraus et al. (2011) [kraus2011]_. 
- **Production**: There are two production function models available: ``neukum`` and ``powerlaw``. The default is ``neukum``, and there are three versions available of this model: ``Moon`` [neukum2001]_, ``Mars`` [ivanov2001-mars]_, and ``Projectile`` [ivanov2001-proj]_.  The version will be selected based on the target body. When the target body is either ``Moon`` or ``Mars``, then one of the two crater-based production functions are selected, and projectile sizes are determined by the Scaling model. For other bodies, ``Projectile`` is chosen, and the crater size is determined by the Scaling model.  
- **Projectile**: There are two projectile velocity and density models available: ``asteroids`` and ``comets``. The default projectile model is determined by the target body. If the body is an inner solar system body, then ``asteroids`` is chosen, otherwise ``comets`` is chosen. The asteroid velocity distributions for the Moon are from Yue et al. (2013) [yue2013]_, and for all other inner planets from Minton et al. (2010) [minton2010]_. Comet velocities are from Zahnle et al. (2003) [zahnle2003]_.
- **Grid**: There are three grid models available: ``icosphere``, ``arbitrary_resolution``, and ``hireslocal``. The default is ``icosphere``, which builds fast an efficient representation of a sphere. The *resolution* of the grid (the number of faces of the mesh) is determined by the formula :math:`20 \times 4^n`, where :math:`n` is given by the argument ``gridlevel`` with a default value of 8.
- **Morphology**: Currently one morphology model is available: ``simplemoon``. This is a model that similar to that used by CTEM. Most of the parameters are taken from Pike (1977) [pike1977]_, except for simple crater profiles, which use a model from Fassett and Thomson (2014) [fassett2014]_. Ejecta blanket scaling is from McGetchin et al. (1973) [mcgetchin1973]_.  

Because the ``Simulation`` class contains all other components, the defaults for all of the components can be viewed from printing a simulation object to the console, or by inspecting the ``cratermaker.yaml`` configuration file of a simulation with no arguments passed to it.

.. ipython:: python
    :okwarning:
    
    import cratermaker as cm
    sim = cm.Simulation()
    print(sim)
