from pathlib import Path
import numpy as np
from numpy.random import Generator
from dataclasses import dataclass
from cratermaker.utils.general_utils import parameter
from typing import Any

@dataclass
class CommonArgs:
    simdir: Path
    rng: Generator
    seed: int | None

class CratermakerBase:
    def __init__(self, 
                 simdir: str | Path = Path.cwd(),
                 rng: Generator | None = None, 
                 seed: int | None = None,
                 **kwargs):
        """
        Initialize the CratermakerBase class.

        Parameters
        ----------
        simdir : str | Path
            The main project simulation directory.
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the seed if it is provided.
        seed : int | None
            The random seed for the simulation if rng is not provided. If None, a random seed is used.
        kwargs : Any
            Additional keyword arguments for subclasses.
        """
        object.__setattr__(self, "_user_defined", set())
        self.simdir = simdir
        self.seed = seed  

        if rng is not None:
            self.rng = rng  # Ignore the seed and use the provided RNG
        else:
            self.rng = np.random.default_rng(seed=seed)

        super().__init__()

    def to_config(self, **kwargs: Any) -> dict[str, Any]:
        def _convert_for_yaml(obj):
            """
            Converts values to types that can be used in yaml.safe_dump. This will convert various types into a format that can be saved in a human-readable YAML file. Therefore, it will ignore anything that cannot be converted into a str, int, float, or bool.
            """
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

        config = _convert_for_yaml({name: getattr(self, name) for name in self._user_defined if hasattr(self, name)})
        return {key: value for key, value in config.items() if value is not None} 

    @parameter
    def simdir(self):
        """
        The main project simulation directory.
        """
        return self._simdir

    @simdir.setter
    def simdir(self, value):
        if value is None:
            self._simdir = Path.cwd()
        else:
            try:
                p = Path(value)
                if not p.is_absolute():
                    p = Path.cwd() / p
                self._simdir = p
                self._simdir.mkdir(parents=True, exist_ok=True)
            except TypeError:
                raise TypeError("simdir must be a path-like object (str, Path, or os.PathLike)")
        
    @parameter
    def seed(self):
        """
        The random seed for the simulation.
        """
        return self._seed

    @seed.setter
    def seed(self, value):
        if value is not None:
            if not isinstance(value, int) or np.isnan(value) or np.isinf(value) or value < 0:
                raise TypeError("seed must be a positive integer")
            self._seed = int(value)
        else:
            self._seed = None

    @property
    def rng(self):
        return self._rng

    @rng.setter
    def rng(self, value):
        if not isinstance(value, Generator):
            raise TypeError("Expected a numpy.random.Generator")
        self._rng = value  

    @property
    def common_args(self) -> CommonArgs:
        return CommonArgs(simdir=self.simdir, rng=self.rng, seed=self.seed)