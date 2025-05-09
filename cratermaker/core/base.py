from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.random import BitGenerator, Generator, RandomState, SeedSequence
from numpy.typing import ArrayLike

from cratermaker.utils.general_utils import parameter


@dataclass
class CommonArgs:
    simdir: Path
    rng: Generator | None
    rng_seed: int | None
    rng_state: dict | None


class CratermakerBase:
    def __init__(
        self,
        simdir: str | Path | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs,
    ):
        """
        Initialize the CratermakerBase class.

        Parameters
        ----------
        simdir : str | Path
            The main project simulation directory. Defaults to the current working directory if None.
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            The rng_rng_seed for the RNG. If None, a new RNG is created.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        **kwargs : Any
            Additional keyword arguments.
        """
        object.__setattr__(self, "_user_defined", set())
        object.__setattr__(self, "_rng", None)
        object.__setattr__(self, "_rng_seed", None)
        object.__setattr__(self, "_rng_state", None)
        object.__setattr__(self, "_simdir", None)
        self.simdir = simdir

        self._rng_seed = rng_seed
        self.rng, self.rng_state = _rng_init(
            rng=rng, rng_seed=rng_seed, rng_state=rng_state
        )

        super().__init__()

    def to_config(
        self, remove_common_args: bool = False, **kwargs: Any
    ) -> dict[str, Any]:
        """
        Converts values to types that can be used in yaml.safe_dump. This will convert various types into a format that can be saved in a human-readable YAML file.

        Parameters
        ----------
        obj : Any
            The object whose attributes will be stored.  It must have a _user_defined attribute.
        remove_common_args : bool, optional
            If True, remove the set of common arguments that are shared among all components of the project from the configuration. Defaults to False.
        **kwargs : Any
            Additional keyword arguments for subclasses.

        Returns
        -------
        dict[str, Any]
            A dictionary of the object's attributes that can be serialized to YAML.
        Notes
        -----
        - The function will ignore any attributes that are not serializable to human-readable YAML. Therefore, it will ignore anything that cannot be converted into a str, int, float, or bool.
        - The function will convert Numpy types to their native Python types.
        """
        return _to_config(self, remove_common_args=remove_common_args, **kwargs)

    @parameter
    def simdir(self):
        """
        The main project simulation directory.

        Returns
        -------
        Path
            The initialized simulation directory as a Path object. Will be a relative path if possible, otherwise will be absolute. If it doesn't exist, it will be created.
        """
        return self._simdir

    @simdir.setter
    def simdir(self, value):
        self._simdir = _simdir_init(value)

    @parameter
    def rng_seed(self):
        """
        The random rng_seed for the simulation RNG.

        Returns
        -------
        int or None
            The integer rng_seed used to initialize the RNG, or None if not set.
        """
        return self._rng_seed

    @rng_seed.setter
    def rng_seed(self, value):
        if value is not None:
            if (
                not isinstance(value, int)
                or np.isnan(value)
                or np.isinf(value)
                or value < 0
            ):
                raise TypeError("rng_seed must be a positive integer")
            self._rng_seed = int(value)
        else:
            self._rng_seed = None

    @property
    def rng(self):
        """
        The random number generator used for stochastic elements of the simulation.

        Returns
        -------
        numpy.random.Generator or None
            The RNG instance, or None if not initialized.
        """
        return self._rng

    @rng.setter
    def rng(self, value):
        self._rng, _ = _rng_init(
            rng=value, rng_seed=self.rng_seed, rng_state=self.rng_state
        )

    @parameter
    def rng_state(self):
        """
        The state of the random number generator.

        Returns
        -------
        dict or None
            A dictionary representing the RNG state, or None if the RNG is not initialized.
        """
        return self.rng.bit_generator.state if self.rng is not None else None

    @rng_state.setter
    def rng_state(self, value):
        _, self._rng_state = _rng_init(
            rng=self.rng, rng_seed=self.rng_seed, rng_state=value
        )

    @property
    def common_args(self) -> CommonArgs:
        return CommonArgs(
            simdir=self.simdir,
            rng=self.rng,
            rng_seed=self.rng_seed,
            rng_state=self.rng_state,
        )


def _rng_init(
    rng: Generator | None = None,
    rng_seed: int
    | ArrayLike
    | SeedSequence
    | BitGenerator
    | Generator
    | RandomState
    | None = None,
    rng_state: dict | None = None,
    **kwargs: Any,
) -> tuple[Generator, dict]:
    """
    Initialize the random number generator (RNG) based on the provided rng_seed.

    Parameters
    ----------
    rng : Generator, optional
        The random number generator to be initialized.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
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
    - If neither rng nor rng_seed is provided, a new RNG will be created using the default rng_seed and the state will be returned.
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


def _simdir_init(simdir: str | Path | None = None, **kwargs: Any) -> Path:
    """
    Initialize the simulation directory.

    Parameters
    ----------
    simdir : str | Path | None
        The main project simulation directory. Defaults to the current working directory if None.

    Returns
    -------
    Path
        The initialized simulation directory as a Path object. Will be a relative path if possible, otherwise will be absolute.
    """

    if simdir is None:
        p = Path.cwd()
    else:
        try:
            p = Path(simdir)
            if not p.is_absolute():
                p = Path.cwd() / p
            p.mkdir(parents=True, exist_ok=True)
            p = p.resolve()
        except TypeError:
            raise TypeError("simdir must be a path-like object (str, Path, or None)")
    try:
        simdir = p.relative_to(Path.cwd())
    except ValueError:
        simdir = p
    return simdir


def _convert_for_yaml(obj):
    if isinstance(obj, dict):
        return {k: _convert_for_yaml(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_convert_for_yaml(v) for v in obj]
    elif isinstance(obj, tuple):
        return tuple(_convert_for_yaml(v) for v in obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    elif isinstance(obj, Path):
        return str(obj)
    elif isinstance(obj, (str, int, float, bool)):
        return obj
    elif obj is None:
        return None
    else:
        return str(obj)


def _to_config(obj, remove_common_args: bool = False, **kwargs: Any) -> dict[str, Any]:
    config = _convert_for_yaml(
        {name: getattr(obj, name) for name in obj._user_defined if hasattr(obj, name)}
    )
    if remove_common_args:
        config = {
            key: value
            for key, value in config.items()
            if key not in obj.common_args.__dict__
        }
    return {key: value for key, value in config.items() if value is not None}
