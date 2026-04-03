r"""
Run a simulation of the Moon with basins emplaced using QuasiMC mode
====================================================================

.. rubric:: By Austin Blevins and David Minton

This example shows how to run a lunar simulation in 'Quasi-Monte Carlo' mode, where the largest craters are read from a csv file that gives their diameters and locations, along with one or more columns that indicate an emplacement time or time range.  For details, see :ref:`ug-production-quasimc` for details.

We will use the |sim.run| method in a Cratermaker Simulation object for an age of 4.31 billion years using the default Neukum production function [#]_. The largest craters (aka basins) will not follow this production function, but instead will be emplaced according to the csv file with given values of N(D), which represents the time needed to produce N craters larger than D diameter in kilometers, per 10⁶ km² of surface area, its age, and/or a stratigraphic sequence number.  Like example 1.2, we also reduce the gridlevel to 6 to speed up the simulation for this example. Even at low resolution, it will take several minutes to run the simulation. We will plot using a reversed grayscale colormap to somewhat mimic the appearence of the lunar surface. The constriaints on basin are provided below, and are similar to those used in Blevins et al. (2025) [#]_, but with a few updates.

Sources
-------

- The oldest crater in the file is South Pole Aitken basin, which we have given an N(20) of 999.8 per 10⁶ km² and a sequence number of 0. This puts it just after our chosen simulation start time in terms of NPF model years.
- The 74 largest and oldest crater names, sizes, and impact locations are those determined using GRAIL-derived gravity data from Neumann et al. (2015) [#]_
- Sequence numbers are determined by stratigraphic relationships from Wilhelms et al. (1987) [#]_ and Fassett et al. (2012) [#]_
- Fecunditatis and Australe North are assigned N(90) values from Evans et al. (2018) [#]_.
- We do not consider the evidence for Procellarum being an impact basin as particularly strong, and so we do not include it.
- Basins with on N(20) values or those with stratigraphic sequences but no ages are taken from Orgel et al. (2018) [#]_.
- Imbrium is set to 3922±12 My based on U-Pb ages of Apollo 14 and 15 impact breccias from Nemchin et al. (2021) [#]_.
- The age of Copernicus of 800±15 My is from Bogard et al. (1994) [#]_, and is based on :sup:`39`\ Ar - :sup:`40`\ Ar ages of Apollo 12 soil samples that are thought to have been reset by ejecta from Copernicus.
- The age of Tycho of 109±4 is based on cosmic ray exposure dating of the South Massif landslide deposit in Taurus Littrow Valley from Apollo 17 by Drozd et al. (1977) [#]_.
- Serenetatis is given an age of 4250 My based on modeling the origin of the Apollo 17 troctolite 76535 (which we think is The Most Interesting Rock from the Moon) by Bjonnes et al. (2015) [#]_. However, this is still a proposed age and because it would represent a date that is not part of the calibration used to develop the Neukum chronology function, then it may not be consistent with its position in the lunar chronology. Alternatively, Orgel et al.'s N(20)=334±73 could be used instead.



References
----------

.. [#] Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55-86. `doi:10.1023/A:1011989004263 <https://doi.org/10.1023/A:1011989004263>`_
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




"""

import cratermaker as cm

simdir = "simdata-1_8"

# Note, that for these examples we pass ask_overwrite=False and reset=True to the Simulation constructor. This will suppress
# prompts that ask the user if they want to overwrite existing files, which would would prevent these examples from running on their
# own when building the documentation pages. Alternatively, calling cm.cleanup(simdir) will remove all pre-existing output files.


# Initialize a quick Moon simulation. We will reduce the resolution to gridlevel 6 and turn off counting to speed up the simulation for this example.
sim = cm.Simulation(
    target="Moon", gridlevel=6, do_counting=False, simdir=simdir, ask_overwrite=False, reset=True, quasimc_file="qmc_input.csv"
)

# When assigning quasimc craters, the size of the smallest Quasi-Monte Carlo crater will be used to set the largest Monte Carlo-generated crater. In our example test case, we have two small young craters (Tycho and Copernicus), so we need to ensure we don't artificially cut off our production population at the size of Tycho.
sim.largest_crater = 200e3
sim.run(age=4310)
sim.show3d(variable_name="face_elevation", cmap="Greys_r")


###############################################################################
# Quasi Monte Carlo input values
# ------------------------------
#
#  The CSV file used in this example is :download:`qmc_input.csv </_static/qmc_input.csv>`
#
# .. csv-table:: qmc_input.csv
#   :file: /_static/qmc_input.csv
#   :widths: auto
#
