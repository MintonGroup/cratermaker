from pathlib import Path
import numpy as np
from numpy.random import Generator
from dataclasses import dataclass
from typing import Any
from cratermaker.utils.general_utils import parameter
from cratermaker.utils.init_helpers import _rng_init

@dataclass
class CommonArgs:
    simdir: Path
    rng: Generator | None
    rng_seed: int | None
    rng_state: dict | None

class CratermakerBase:
    def __init__(self, 
                 simdir: str | Path = Path.cwd(),
                 rng: Generator | None = None, 
                 rng_seed: int | None = None,
                 rng_state: dict | None = None,
                 **kwargs):
        """
        Initialize the CratermakerBase class.

        Parameters
        ----------
        simdir : str | Path
            The main project simulation directory.
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
        rng_seed : Any type allowed by the seed argument of numpy.random.Generator, optional
            The rng_rng_seed for the RNG. If None, a new RNG is created.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        kwargs : Any
            Additional keyword arguments for subclasses.
        """
        object.__setattr__(self, "_user_defined", set())
        self.simdir = simdir

        self._rng_seed = rng_seed
        self.rng, self.rng_state  = _rng_init(rng=rng, rng_seed=rng_seed, rng_state=rng_state)

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
    def rng_seed(self):
        """
        The random rng_seed for the simulation.
        """
        return self._rng_seed

    @rng_seed.setter
    def rng_seed(self, value):
        if value is not None:
            if not isinstance(value, int) or np.isnan(value) or np.isinf(value) or value < 0:
                raise TypeError("rng_seed must be a positive integer")
            self._rng_seed = int(value)
        else:
            self._rng_seed = None

    @property
    def rng(self):
        return self._rng

    @rng.setter
    def rng(self, value):
        if not isinstance(value, Generator):
            raise TypeError("Expected a numpy.random.Generator")
        self._rng = value  

    @parameter 
    def rng_state(self):
        """
        The state of the random number generator.
        """
        return self._rng_state
    
    @rng_state.setter
    def rng_state(self, value):
        if value is not None:
            try:
                tmp = np.random.default_rng()
                tmp.bit_generator.state = value
            except Exception as e:
                raise TypeError("rng_state is not a valid bit_generator state") from e
            self._rng_state = value
        else:
            self._rng_state = None

    @property
    def common_args(self) -> CommonArgs:
        return CommonArgs(simdir=self.simdir, rng=self.rng, rng_seed=self.rng_seed, state=self..state)