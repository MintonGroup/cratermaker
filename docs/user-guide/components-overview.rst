.. currentmodule:: cratermaker

.. _ug-components-overview:

#########################
Components of Cratermaker
#########################

Each component of Cratermaker is built around a set of classes, which include:

- **Target**: Models the properties of the target body, including its size, material properties, and surface gravity. 
- **Projectile**: Models the properties of projectiles, including their size, velocity, and density. 
- **Scaling**: Contains the model for relating projectile size and crater size. 
- **Production**: Manages the production function, which determines the number of craters/projectiles produced over time and surface area. 
- **Surface**: Represents the surface of the target body, including its resolution, topography, material properties. 
- **Morphology**: Used to construct the three-dimensional representation of a crater onto the surface. 
- **Crater**: Represents a single crater and the properties of the projectile that formed it. 
- **Simulation**: This is main class that can be used to run a comprehensive cratered landscape evolution model. 

The class definition of each component contains a special construction method called ``.maker()``, which is used to create the component. This method can be used to create a component without having to create the entire simulation. While some components require others to function (for instance, a Scaling model needs both a Target and Projectile component in order to perform the scaling computations), when creating such an object, any missing components are built internally using resonable defaults (see :ref:`ug-defaults` for details).

.. _ug-target:
.. include:: target.rst
.. _ug-projectile:
.. include:: projectile.rst
.. _ug-scaling:
.. include:: scaling.rst
.. _ug-production:
.. include:: production.rst
.. _ug-surface:
.. include:: surface.rst
.. _ug-morphology:
.. include:: morphology.rst
.. _ug-crater:
.. include:: crater.rst
.. _ug-defaults:
.. include:: defaults.rst

.. _ug-references:

References
==========

.. [richardson2009] Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029
.. [holsapple1993] Holsapple, K.A., 1993. The scaling of impact processes in planetary sciences 21, 333-373. https://doi.org/10.1146/annurev.ea.21.050193.002001
.. [kraus2011] Kraus, R.G., Senft, L.E., Stewart, S.T., 2011. Impacts onto H2O ice: Scaling laws for melting, vaporization, excavation, and final crater size. Icarus 214, 724-738. https://doi.org/10.1016/j.icarus.2011.05.016
.. [neukum2001] Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55–86. https://doi.org/10.1023/A:1011989004263
.. [ivanov2001-proj] Ivanov, B.A., Neukum, G., Wagner, R., 2001. Size-Frequency Distributions of Planetary Impact Craters and Asteroids, in: Collisional Processes in the Solar System. Springer Netherlands, Dordrecht, pp. 1–34. https://doi.org/10.1007/978-94-010-0712-2_1
.. [ivanov2001-mars] Ivanov, B.A., 2001. Mars/Moon Cratering Rate Ratio Estimates. *Space Science Reviews*, 96, 87-104. https://doi.org/10.1023/A:1011941121102
.. [yue2013] Yue, Z., Johnson, B.C., Minton, D.A., Melosh, H.J., Di, K., Hu, W., Liu, Y., 2013. Projectile remnants in central peaks of lunar impact craters. Nature Geosci 6, 435 EP-. https://doi.org/10.1038/ngeo1828
.. [minton2010] Minton, D.A., Malhotra, R., 2010. Dynamical erosion of the asteroid belt and implications for large impacts in the inner Solar System. Icarus 207, 744–757. https://doi.org/10.1016/j.icarus.2009.12.008
.. [zahnle2003] Zahnle, K., Schenk, P., Levison, H., Dones, L., 2003. Cratering rates in the outer Solar System. Icarus 163, 263–289. https://doi.org/10.1016/S0019-1035(03)00048-4
.. [pike1977] Pike, R.J., 1977. Size-dependence in the shape of fresh impact craters on the moon. Presented at the In: Impact and explosion cratering: Planetary and terrestrial implications; Proceedings of the Symposium on Planetary Cratering Mechanics, pp. 489–509.
.. [fassett2014] Fassett, C.I., Thomson, B.J., 2014. Crater degradation on the lunar maria: Topographic diffusion and the rate of erosion on the Moon. J. Geophys. Res. 119, 2014JE004698-2271. https://doi.org/10.1002/2014JE004698
.. [mcgetchin1973] McGetchin, T.R., Settle, M., Head, J.W., 1973. Radial thickness variation in impact crater ejecta: Implications for lunar basin deposits. Earth Planet. Sci. Lett. 20, 226–236. https://doi.org/10.1016/0012-821X(73)90162-3
