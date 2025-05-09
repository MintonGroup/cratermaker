.. _ug-defaults:

Default Behavior
================
As Cratermaker is designed to be easy to use, all of its component classes are built to be invoked without any arguments. When a ``Simulation`` object is created as in the examples above, a set of component objects are created with their default parameters. It is important to understand what the default behavior of these classes is, as this will affect the results of your simulation.

- **Target**: The default target is the Moon. There are a number of known bodies that can be selected as targets, including some jovian satellites, and small bodies, but currently few of them have been tested. 
- **Scaling**: The default scaling model is called ``default``, as it is similar to the one used by the Cratered Terrain Evolution Model (CTEM) that is a progenitor to Cratermaker [#]_. The projectile to transient scaling model is mostly based on Holsapple (1993) [#]_ with some additional scaling parameters for ice given by Kraus et al. (2011) [#]_. 
- **Production**: There are two production function models available: ``neukum`` and ``powerlaw``. The default is ``neukum``, and there are three versions available of this model: ``Moon`` [#]_, ``Mars`` [#]_, and ``Projectile`` [#]_.  The version will be selected based on the target body. When the target body is either ``Moon`` or ``Mars``, then one of the two crater-based production functions are selected, and projectile sizes are determined by the Scaling model. For other bodies, ``Projectile`` is chosen, and the crater size is determined by the Scaling model.  
- **Projectile**: There are two projectile velocity and density models available: ``asteroids`` and ``comets``. The default projectile model is determined by the target body. If the body is an inner solar system body, then ``asteroids`` is chosen, otherwise ``comets`` is chosen. The asteroid velocity distributions for the Moon are from Yue et al. (2013) [#]_, and for all other inner planets from Minton et al. (2010) [#]_. Comet velocities are from Zahnle et al. (2003) [#]_.
- **Grid**: There are three grid models available: ``icosphere``, ``arbitrary_resolution``, and ``hireslocal``. The default is ``icosphere``, which builds fast an efficient representation of a sphere. The *resolution* of the grid (the number of faces of the mesh) is determined by the formula :math:`20 \times 4^n`, where :math:`n` is given by the argument ``gridlevel`` with a default value of 8.
- **Morphology**: Currently one morphology model is available: ``simplemoon``. This is a model that similar to that used by CTEM. Most of the parameters are taken from Pike (1977) [#]_, except for simple crater profiles, which use a model from Fassett and Thomson (2014) [#]_. Ejecta blanket scaling is from McGetchin et al. (1973) [#]_.  

Because the ``Simulation`` class contains all other components, the defaults for all of the components can be viewed from printing a simulation object to the console, or by inspecting the ``cratermaker.yaml`` configuration file of a simulation with no arguments passed to it.

.. ipython:: python
    :okwarning:
    
    import cratermaker as cm
    sim = cm.Simulation()
    print(sim)


References
----------

.. [#] Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029
.. [#] Holsapple, K.A., 1993. The scaling of impact processes in planetary sciences 21, 333-373. https://doi.org/10.1146/annurev.ea.21.050193.002001
.. [#] Kraus, R.G., Senft, L.E., Stewart, S.T., 2011. Impacts onto H2O ice: Scaling laws for melting, vaporization, excavation, and final crater size. Icarus 214, 724-738. https://doi.org/10.1016/j.icarus.2011.05.016
.. [#] Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55-86. https://doi.org/10.1023/A:1011989004263
.. [#] Ivanov, B.A., 2001. Mars/Moon Cratering Rate Ratio Estimates. *Space Science Reviews*, 96, 87-104. https://doi.org/10.1023/A:1011941121102
.. [#] Ivanov, B.A., Neukum, G., Wagner, R., 2001. Size-Frequency Distributions of Planetary Impact Craters and Asteroids, in: Collisional Processes in the Solar System. Springer Netherlands, Dordrecht, pp. 1-34. https://doi.org/10.1007/978-94-010-0712-2_1
.. [#] Yue, Z., Johnson, B.C., Minton, D.A., Melosh, H.J., Di, K., Hu, W., Liu, Y., 2013. Projectile remnants in central peaks of lunar impact craters. Nature Geosci 6, 435 EP-. https://doi.org/10.1038/ngeo1828
.. [#] Minton, D.A., Malhotra, R., 2010. Dynamical erosion of the asteroid belt and implications for large impacts in the inner Solar System. Icarus 207, 744-757. https://doi.org/10.1016/j.icarus.2009.12.008
.. [#] Zahnle, K., Schenk, P., Levison, H., Dones, L., 2003. Cratering rates in the outer Solar System. Icarus 163, 263-289. https://doi.org/10.1016/S0019-1035(03)00048-4
.. [#] Pike, R.J., 1977. Size-dependence in the shape of fresh impact craters on the moon. Presented at the In: Impact and explosion cratering: Planetary and terrestrial implications; Proceedings of the Symposium on Planetary Cratering Mechanics, pp. 489-509.
.. [#] Fassett, C.I., Thomson, B.J., 2014. Crater degradation on the lunar maria: Topographic diffusion and the rate of erosion on the Moon. J. Geophys. Res. 119, 2014JE004698-2271. https://doi.org/10.1002/2014JE004698
.. [#] McGetchin, T.R., Settle, M., Head, J.W., 1973. Radial thickness variation in impact crater ejecta: Implications for lunar basin deposits. Earth Planet. Sci. Lett. 20, 226-236. https://doi.org/10.1016/0012-821X(73)90162-3
