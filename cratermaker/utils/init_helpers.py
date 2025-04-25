import numpy as np
from numpy.random import Generator, SeedSequence, BitGenerator, RandomState
from numpy.typing import ArrayLike
from pathlib import Path
from cratermaker import Target, Surface, Crater
from cratermaker import ScalingModel, get_scaling_model, available_scaling_models
from cratermaker import ProductionModel, get_production_model, available_production_models
from cratermaker import MorphologyModel, get_morphology_model, available_morphology_models
from cratermaker import ImpactorModel, get_impactor_model, available_impactor_models
from cratermaker import GridMaker, get_grid_type, available_grid_types


def _rng_init(rng: Generator | None = None, 
              rng_seed:  int | ArrayLike | SeedSequence | BitGenerator | Generator | RandomState | None = None,
              rng_state: dict | None = None,
              **kwargs) -> tuple[Generator, dict]:
    """
    Initialize the random number generator (RNG) based on the provided rng_seed.

    Parameters
    ----------
    rng : Generator, optional
        The random number generator to be initialized.
    rng_seed : Any type allowed by the seed argument of numpy.random.Generator, optional
        The rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        Set the rng_state of the RNG. 

    Returns
    -------
    Generator
        The initialized RNG instance.
    dict
        The state of the RNG.

    Notes
    -----
    - If rng is provided, it will be used directly and any input rng_seed or rng_state will be ignored.
    - If both rng_seed and rng_state are provided, rng_state will take precedence.
    - If neither rng nor rng_seed is provided, a new RNG will be created using the default seed and the state will be returned.
    """
    if rng is not None:
        if not isinstance(rng, Generator):
            raise ValueError("rng must be a numpy.random.Generator instance.")
    else:
        if rng_seed is None: 
            rng = np.random.default_rng()
        if rng_state is not None:
            if not isinstance(rng_state, dict):
                raise ValueError("rng_state must be a dictionary.")
            try:
                rng.bit_generator.state = rng_state
            except Exception as e:
                raise ValueError("Invalid rng_state provided.") from e
        elif rng_seed is not None:
            try:
                rng = np.random.default_rng(seed=rng_seed)
            except Exception as e:
                raise ValueError("Invalid rng_seed provided.") from e

    rng_state = rng.bit_generator.state
    return rng, rng_state
    


