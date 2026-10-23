"""
Perform tasks on running simulations with save actions
======================================================

.. rubric:: By David Minton

When running simulations, sometimes you want to have certain tasks performed on the intermediate results in order to track progress or perform analysis concurrently with the Simulation. This can be accomplished using a property called `save_actions`.  This example shows how you can write a custom function to perform an analysis task on a running simulation and call it using save actions.

"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cratermaker import Simulation

simdir = "simdata-1_10"

sim = Simulation(
    surface="hireslocal",
    local_location=(0, 0),
    pix=400.0,
    local_radius=100.0e3,
    simdir=simdir,
    ask_overwrite=False,
    reset=True,
    rng_seed=87398452,
)


def plot_sfd(interval, **kwargs):
    if interval == 0 or sim.n_observed == 0:
        return
    fig, ax = plt.subplots(figsize=(3, 5), layout="constrained")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("$\\mathregular{N_{>D}}$")
    ax.set_xlabel("Diameter (km)")
    ax.set_title(f"Time: {sim.time:.0f} My bp")
    ax.set_xlim(1, 20.0)
    ax.set_ylim(0.1, 40)
    ds = sim.counting.to_xarray(sim.observed)
    ds = ds.sortby(-ds.diameter)
    n = np.arange(1, len(ds.diameter) + 1)
    nsig = np.sqrt(n)
    ax.errorbar(ds.diameter * 1e-3, n, yerr=nsig, label="Observed")
    ax.legend(loc="lower left")
    # plt.savefig(f"sfd{interval:02d}.pdf")
    plt.show()
    return


sim.add_save_action({plot_sfd: {}})

sim.run(age=3700, ninterval=5)
