"""
Run a simulation of the Moon with basins emplaced using QuasiMC mode
====================================================================

.. rubric:: By Austin Blevins and David Minton

This example shows how to run a lunar simulation in 'Quasi-Monte Carlo' mode, where the largest craters are read from a csv file that gives their diameters and locations, along with one or more columns that indicate an emplacement time or time range.  For details, see :ref:`ug-production-quasimc` for details.

The emplacement time of a Quasi-Monte Carlo crater is determined by some combination of the following columns.

production_D
    The `D` value (in km) used in the N(D) convention (e.g. 1 would be used to indicate N(1), 20 for N(20), etc.)

production_N
    The N(D) value in units of craters per 10⁶ km².

production_N_stdev
    The 1-sigma standard deviation of a range of N(D) values, given in units of craters per 10⁶ km²

production_time
    The mean emplacement time (in My before present)

production_time_stdev
    The 1-sigma standard deviation of possible time values (in My before present)

production_sequence
    An integer value that indicates the emplacement sequence relative to other craters with a production_sequence value. For example, a crater with production_sequence=5 would be emplaced after a crater with production_sequence=4, but before a crater with production_sequence=6. Craters can have a production_sequence value even if they don't have a production_time or production_N value, as long as at least 1 crater with the lowest production_sequence value has either a production_time or production_N value. Otherwise, the production_N is interpolated based on the production_N/_time values of the neighboring sequences.

If both N(D) and time values are provided, the N(D) values take precedence, and if stdev values are ommitted, this is interpreted as 0.

We will use the |sim.run| method in a Cratermaker Simulation object for an age of 4.31 billion years using the default Neukum production function [#]_. The largest craters (aka basins) will not follow this production function, but instead will be emplaced according to the csv file. The age of South Pole Aitken basin is set to simulation start, and Imbrium is set to 3.9 billion years ago. The other basins can vary based on N(20) values measured in Orgel et al. (2018) [#]_. The 74 basins catalogued by Neumann et al. (2015) [#]_ are included in this file.

Like example 1.2, we also reduce the gridlevel to 6 to speed up the simulation for this example. Even at low resolution, it will take several minutes to run the simulation. We will plot using a reversed grayscale colormap to somewhat mimic the appearence of the lunar surface.


References
----------

.. [#] Neukum, G., Ivanov, B.A., Hartmann, W.K., (2001) Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55-86. `doi: 10.1023/A:1011989004263 <https://doi.org/10.1023/A:1011989004263>`_
.. [#] Orgel, C., Michael, G., Fassett, C.I., van der Bogert, C.H., Riedel, C., Kneissl, T., Hiesinger, H., 2018. Ancient Bombardment of the Inner Solar System: Reinvestigation of the “Fingerprints” of Different Impactor Populations on the Lunar Surface. J. Geophys. Res. 123, 748–762. `doi: 10.1002/2017JE005451 <https://doi.org/10.1002/2017JE005451>`_
.. [#] Neumann, G.A., Zuber, M.T., Wieczorek, M.A., Head, J.W., Baker, D.M.H., Solomon, S.C., Smith, D.E., Lemoine, F.G., Mazarico, E., Sabaka, T.J., Goossens, S.J., Melosh, H.J., Phillips, R.J., Asmar, S.W., Konopliv, A.S., Williams, J.G., Sori, M.M., Soderblom, J.M., Miljković, K., Andrews-Hanna, J.C., Nimmo, F., Kiefer, W.S., 2015. Lunar impact basins revealed by Gravity Recovery and Interior Laboratory measurements. Science Advances 1, e1500852–e1500852. `doi: 10.1126/sciadv.1500852 <https://doi.org/10.1126/sciadv.1500852>`_
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

# When assigning quasimc craters, the size of the smallest Quasi-Monte Carlo crater will be used to set the largest Monte Carlo-generated crater. In our example test case, we have two small young craters (Tycho and Copernicus), so we need to ensure we don't artificially cut off our production population at the size of Tycho.
sim.largest_crater = 200e3
sim.run(age=4310)
sim.show3d(variable_name="face_elevation", cmap="Greys_r")
