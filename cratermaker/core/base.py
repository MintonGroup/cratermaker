from __future__ import annotations

import importlib
import pkgutil
import shutil
import sys
from abc import ABC
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
    """
    Base class for the Cratermaker project.

    Parameters
    ----------
    simdir : str | Path
        The main project simulation directory. Default is the current working directory if None.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments.
    """

    def __init__(
        self,
        simdir: str | Path | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs,
    ):
        object.__setattr__(self, "_user_defined", set())
        object.__setattr__(self, "_rng", None)
        object.__setattr__(self, "_rng_seed", None)
        object.__setattr__(self, "_rng_state", None)
        object.__setattr__(self, "_simdir", None)
        object.__setattr__(self, "_output_dir_name", None)
        object.__setattr__(self, "_output_file_pattern", [])

        self.simdir = simdir

        self._rng_seed = rng_seed
        self.rng, self.rng_state = _rng_init(rng=rng, rng_seed=rng_seed, rng_state=rng_state)

        super().__init__()

    def to_config(self, remove_common_args: bool = False, **kwargs: Any) -> dict[str, Any]:
        """
        Converts values to types that can be used in yaml.safe_dump. This will convert various types into a format that can be saved in a human-readable YAML file.

        Parameters
        ----------
        obj : Any
            The object whose attributes will be stored.  It must have a _user_defined attribute.
        remove_common_args : bool, optional
            If True, remove the set of common arguments that are shared among all components of the project from the configuration. Default is False.
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

    def reprocess(self, interval_range: ArrayLike = (-1), **kwargs: Any) -> None:
        """
        Reprocess the output from this component.

        This is a placeholder method that can be overridden by subclasses to implement specific reprocessing logic.

        Parameters
        ----------
        interval_range : ArrayLike, optional
            The range of intervals to reprocess. Default is (-1), which indicates the last interval.
        """
        pass

    def has_output(self, **kwargs: Any) -> bool:
        """
        Check if the component has any output files in its output directory.

        Returns
        -------
        list[Path]
            A list of Path objects representing the files that would be removed during a reset operation. Returns an empty list if no files found
        """
        if self.output_dir is None or not self.output_file_pattern:
            return []
        if not self.output_dir.exists():
            return []
        files_to_remove = []
        for pattern in self.output_file_pattern:
            files_to_remove.extend(self.output_dir.glob(pattern))
        if files_to_remove:
            return files_to_remove
        else:
            return []

    def reset(self, ask_overwrite: bool = False, files_to_remove: list[Path | str] | None = None, **kwargs: Any) -> None:
        """
        Reset the component by removing its output files.

        Parameters
        ----------
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before deleting files. Default is False.
        files_to_remove : list[Path | str], optional
            If set, this is the list of files that will be removed. If not set, then the removed files will be determined using the component's output_dir and output_file_pattern attributes.
        **kwargs : Any
            Additional keyword arguments for subclasses.
        """
        if files_to_remove is None:
            files_to_remove = self.has_output()
        if files_to_remove:
            if ask_overwrite:
                print(f"The following files will be deleted in {self.output_dir}:")
                for f in files_to_remove:
                    print(f"  {f}")
                print("To disable this message, pass `ask_overwrite=False` to this function.")
                response = input(f"Are you sure you want to delete {len(files_to_remove)} files in {self.output_dir}? [y/N]: ")
                if response.lower() != "y":
                    raise RuntimeError("User aborted the reset operation.")
            for file_path in files_to_remove:
                try:
                    p = Path(file_path)
                    if p.is_dir():
                        shutil.rmtree(p)
                    else:
                        p.unlink()
                except Exception as e:
                    print(f"Error removing file {file_path}: {e}", file=sys.stderr)

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
        if isinstance(value, Path):
            self._simdir = value
        elif isinstance(value, (str | None)):
            self._simdir = _simdir_init(value)

    @property
    def output_dir(self) -> Path | None:
        """
        The output directory for a component. If None, the component does not have an output directory set.
        """
        if self._output_dir_name is None:
            return None
        output_dir = self.simdir / self._output_dir_name
        if not output_dir.exists():
            try:
                output_dir.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                raise RuntimeError(f"Could not create output directory at {output_dir}") from e
        return output_dir

    @property
    def output_file_pattern(self) -> list[str]:
        """
        Return a list of file patterns that this component will generate when saved.

        This is used during the reset operation to remove old files generated by this component if requested.

        Returns
        -------
        list[str]
            A list of file patterns (e.g., ["*.tiff", "*.gpkg"]) that this component will generate when saved.
            If the component does not generate any files, return an empty list.
        """
        return self._output_file_pattern

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
            if not isinstance(value, (int | np.integer)) or np.isnan(value) or np.isinf(value) or value < 0:
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
        self._rng, _ = _rng_init(rng=value, rng_seed=self.rng_seed, rng_state=self.rng_state)

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
        _, self._rng_state = _rng_init(rng=self.rng, rng_seed=self.rng_seed, rng_state=value)

    @property
    def common_args(self) -> CommonArgs:
        return CommonArgs(
            simdir=self.simdir,
            rng=self.rng,
            rng_seed=self.rng_seed,
            rng_state=self.rng_state,
        )


class ComponentBase(CratermakerBase, ABC):
    """
    Base class for components of the Cratermaker project.

    Defines the common parameters and methods for all components in the Cratermaker project, including the maker class that is used to select the correct component from user arguments.

    Parameters
    ----------
    **kwargs : Any
        Additional keyword arguments.
    """

    _registry: dict[str, type[ComponentBase]] = {}

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)

    def __str__(self) -> str:
        # Return just the name of the class
        base_class = type(self).__mro__[1].__name__
        return f"<{base_class}: {self._component_name}>"

    @classmethod
    def maker(
        cls,
        component: str | type[ComponentBase] | ComponentBase | None = None,
        **kwargs: Any,
    ) -> ComponentBase:
        """
        Initialize a component model with the given name or instance.

        Parameters
        ----------
        component : str or ComponentBase or None
            The name of the component to use, or an instance of ComponentBase. If None, it choose a default component.
        kwargs : Any
            Additional keyword arguments to pass to the component model constructor.

        Returns
        -------
        component
            An instance of the specified component model.

        Raises
        ------
        KeyError
            If the specified component model name is not found in the registry.
        TypeError
            If the specified component model is not a string or a subclass of component.
        """
        if component is None:
            component = cls.available()[0]  # Default to the first available component
        if isinstance(component, str):
            if component not in cls.available():
                raise KeyError(f"Unknown component model: {component}. Available models: {cls.available()}")
            return cls._registry[component](**kwargs)
        elif isinstance(component, type) and issubclass(component, ComponentBase):
            return component(**kwargs)
        elif isinstance(component, ComponentBase):
            return component
        else:
            raise TypeError(f"component must be a string or a subclass of component, not {type(component)}")

    @parameter
    def name(self):
        """
        The registered name of this scaling model set by the @ComponentBase.register decorator.
        """
        return self._component_name

    @classmethod
    def register(cls, name: str):
        """
        Class decorator to register a component model component under the given key.
        """

        def decorator(subcls):
            subcls._component_name = name
            subcls._registry[name] = subcls
            return subcls

        return decorator

    @classmethod
    def available(cls) -> list[str]:
        """Return list of all registered catalogue names."""
        return list(cls._registry.keys())


def import_components(package_name: str, package_path: list[str]) -> None:
    """
    Import all modules of a component package, optionally skipping private modules.

    Parameters
    ----------
    package_name : str
        The full name of the package (e.g., "cratermaker.components.component").
    package_path : list[str]
        The __path__ attribute of the package (usually just __path__).
    """
    for _, module_name, _ in pkgutil.iter_modules(package_path):
        importlib.import_module(f"{package_name}.{module_name}")


def _rng_init(
    rng: Generator | None = None,
    rng_seed: int | ArrayLike | SeedSequence | BitGenerator | Generator | RandomState | None = None,
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
        The main project simulation directory. Default is the current working directory if None.

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
        except TypeError as e:
            raise TypeError("simdir must be a path-like object (str, Path, or None)") from e
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
    config = _convert_for_yaml({name: getattr(obj, name) for name in obj._user_defined if hasattr(obj, name)})
    if remove_common_args:
        config = {key: value for key, value in config.items() if key not in obj.common_args.__dict__}
    return {key: value for key, value in config.items() if value is not None}
