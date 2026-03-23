"""
Run a simulation of the Moon with basins emplaced using QuasiMC mode
====================================================================

.. rubric:: By Austin Blevins and David Minton

This example shows how to run a lunar simulation in 'Quasi-Monte Carlo' mode, where the largest craters are read from a csv file that gives their diameters and locations, along with one or more columns that indicate an emplacement time or time range.

The emplacement time of a Quasi-Monte Carlo crater is determined by some combination of the following columns.

production_D
    The `D` value (in km) used in the N(D) convention (e.g. 1 would be used to indicate N(1), 20 for N(20), etc.)

production_N_low
    The low end of a range of N(D) values, given in units of craters per 10⁶ km²

production_N_high
    The high end of a range of N(D) values,  given in units of craters 1 million sq. km.

production_N
    Equivalent to setting the production_N_high and production_N_low to the same value.

production_time_low
    The low end of a range of time values (in My before present)

production_time_high
    The high end of a range of time values (in My before present)

production_time
    Equivalent to setting production_time_low and production_time_high to be the same.

If both N(D) and time values are provided, the N(D) values take precedence, and if ranges are provided with singular values (etc. both production_time and production_time_low/high are included), the range values take precedence over the singular values. The included file `qmc_input.csv` contains inputs in all of the possible formats discussed above.

We will use the ``run`` method in a Cratermaker Simulation object for an age of 4.31 billion years using the default Neukum production function [#]_. The largest craters (aka basins) will not follow this production function, but instead will be emplaced according to the csv file. The age of South Pole Aitken basin is set to simulation start, and Imbrium is set to 3.9 billion years ago. The other basins can vary based on N(20) values measured in Orgel et al. (2018) [#]_. The 74 basins catalogued by Neumann et al. (2015) [#]_ are included in this file.

Like example 1.2, we also reduce the gridlevel to 6 to speed up the simulation for this example. Even at low resolution, it will take several minutes to run the simulation. We will plot using a reversed grayscale colormap to somewhat mimic the appearence of the lunar surface.


References
----------

.. [#] Neukum, G., Ivanov, B.A., Hartmann, W.K., (2001) Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55-86. `doi: 10.1023/A:1011989004263 <https://doi.org/10.1023/A:1011989004263>`
.. [#] Orgel, C., Michael, G., Fassett, C.I., van der Bogert, C.H., Riedel, C., Kneissl, T., Hiesinger, H., 2018. Ancient Bombardment of the Inner Solar System: Reinvestigation of the “Fingerprints” of Different Impactor Populations on the Lunar Surface. J. Geophys. Res. 123, 748–762. `doi: 10.1002/2017JE005451 <https://doi.org/10.1002/2017JE005451>`
.. [#] Neumann, G.A., Zuber, M.T., Wieczorek, M.A., Head, J.W., Baker, D.M.H., Solomon, S.C., Smith, D.E., Lemoine, F.G., Mazarico, E., Sabaka, T.J., Goossens, S.J., Melosh, H.J., Phillips, R.J., Asmar, S.W., Konopliv, A.S., Williams, J.G., Sori, M.M., Soderblom, J.M., Miljković, K., Andrews-Hanna, J.C., Nimmo, F., Kiefer, W.S., 2015. Lunar impact basins revealed by Gravity Recovery and Interior Laboratory measurements. Science Advances 1, e1500852–e1500852. `doi: 10.1126/sciadv.1500852 <https://doi.org/10.1126/sciadv.1500852>`
"""

import cratermaker as cm

simdir = "simdata-1_7"

# Note, that for these examples we pass ask_overwrite=False and reset=True to the Simulation constructor. This will suppress
# prompts that ask the user if they want to overwrite existing files, which would would prevent these examples from running on their
# own when building the documentation pages. Alternatively, calling cm.cleanup(simdir) will remove all pre-existing output files.


# Initialize a quick Moon simulation. We will reduce the resolution to gridlevel 6 and turn off counting to speed up the simulation for this example.
sim = cm.Simulation(
    target="Moon", gridlevel=6, do_counting=False, simdir=simdir, ask_overwrite=False, reset=True, quasimc_file="qmc_input.csv"
)

sim.run(age=4310)
sim.show3d(variable_name="face_elevation", cmap="Greys_r")
