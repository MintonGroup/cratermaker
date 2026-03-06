from __future__ import annotations

import importlib
import pkgutil
import shutil
import sys
from abc import ABC
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from numpy.random import BitGenerator, Generator, RandomState, SeedSequence
from numpy.typing import ArrayLike
from tqdm import tqdm

from cratermaker.utils.general_utils import parameter


@dataclass
class CommonArgs:
    """
    Stores the common arguments that are shared among all components of the project.

    This is used to pass these common arguments to the maker functions of each component when they are initialized from other components. For example, the Scaling component will store an instance of Target and Projectile. If defaults of these components are used, then the Scaling.maker will pass the arguments stored here to the Target.maker and Projectile.maker functions when it initializes the Target and Projectile objects. This way, if the user has specified any non-default arguments for the simulation (e.g. a custom simdir or rng_seed), then those arguments will be passed to the Target and Projectile makers when they are initialized from the Scaling maker, ensuring that all components of the simulation are initialized with the same common arguments.
    """

    simdir: Path
    """The root directory of the simulation."""
    rng: Generator | None
    """The random number generat or object."""
    rng_seed: int | None
    """The random number seed."""
    rng_state: dict | None
    """The current state of the random number generator."""
    ask_overwrite: bool = True
    """Flag indicating whether to prompt the user before overwriting a file."""


class CratermakerBase:
    """
    Base class for the Cratermaker project.

    Parameters
    ----------
    simdir : str | Path
        |simdir|
    rng : numpy.random.Generator | None
        |rng|
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        |rng_seed|
    rng_state : dict, optional
        |rng_state|
    ask_overwrite : bool, optional
        If True, prompt the user for confirmation before overwriting any existing files when saving output. Default is True.
    **kwargs : Any
        |kwargs|
    """

    def __init__(
        self,
        simdir: str | Path | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        ask_overwrite: bool | None = None,
        **kwargs,
    ):
        object.__setattr__(self, "_user_defined", set())
        object.__setattr__(self, "_rng", None)
        object.__setattr__(self, "_rng_seed", None)
        object.__setattr__(self, "_rng_state", None)
        object.__setattr__(self, "_simdir", None)
        object.__setattr__(self, "_output_dir_name", None)
        object.__setattr__(self, "_output_file_pattern", [])
        object.__setattr__(self, "_output_file_prefix", None)
        object.__setattr__(self, "_output_file_extension", "nc")
        object.__setattr__(self, "_output_image_file_extension", "png")
        object.__setattr__(self, "_export_dir_name", "export")
        object.__setattr__(self, "_ask_overwrite", None)
        object.__setattr__(self, "_save_actions", [])
        object.__setattr__(self, "_driver_to_extension_map", {})

        self.simdir = simdir
        self.ask_overwrite = ask_overwrite

        self._rng_seed = rng_seed
        self.rng, self.rng_state = _rng_init(rng=rng, rng_seed=rng_seed, rng_state=rng_state)

        super().__init__()

    def reset(
        self,
        files_to_remove: list[Path | str] | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Reset the component by removing its output files.

        Parameters
        ----------
        files_to_remove : list[Path | str], optional
            If set, this is the list of files that will be removed. If not set, then the removed files will be determined using the component's output_dir and output_file_pattern attributes.
        **kwargs : Any
            |kwargs|
        """
        if files_to_remove is None:
            files_to_remove = self.saved_output_files()
        if files_to_remove and not self._overwrite_check(files_to_remove):
            print("Reset operation cancelled by user.")
            return

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
            |kwargs|

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

    def save(
        self,
        interval: int | None = None,
        filename: str | Path | None = None,
        skip_actions: bool = False,
        **kwargs: Any,
    ) -> None:
        """
        Runs the action methods specified in the save_actions property of this component, which should be a dictionary where the keys are the names of actions that can be performed on this component (e.g. "plot") when calling the save function and the values are the arguments that should be passed to that action when it is called by the save function.

        Parameters
        ----------
        interval : int or None, optional
            The interval number to save. If None, save without an interval number. Default is None.
        filename : str or Path or None, optional
            The base filename to save to, without the output directory or file extension. If None, the filename will be generated using the output_file_prefix and output_file_extension properties. Default is None.
        skip_actions : bool, optional
            If True, skip running the actions specified in the save_actions property and just return without doing anything. Default is False.
        **kwargs : Any
            |kwargs|

        Notes
        -----
        - skip_actions should be set to True if the save method is called from inside another method of the same object that could potentially be included in a save_actions entry. This will prevent recursive calls to save.
        """
        if skip_actions:
            return

        if filename is None:
            filename = self.output_filename(interval=interval)

        for entry in self.save_actions:
            for action, action_kwargs in entry.items():
                if hasattr(self, action):
                    action_method = getattr(self, action)
                    if callable(action_method):
                        args = {**action_kwargs, **kwargs, "interval": interval}
                        action_method(**args)
                    else:
                        raise ValueError(f"{action} is not a valid action for {self.name}")
                else:
                    raise ValueError(f"{action} is not a valid action for {self.name}")
        return

    def export(self, **kwargs: Any) -> None:
        """
        Export the component data to the specified format.

        This is a stub that acts as a pass-through for components that don't have an export method defined.

        Parameters
        ----------
        **kwargs : Any
            |kwargs|
        """
        pass

    def from_file(self, filename: str | Path, **kwargs: Any) -> None:
        """
        Load the component data from a file.

        This is a stub that acts as a pass-through for components that don't have a from_file method defined.

        Parameters
        ----------
        filename : str or Path
            The path to the file to load the component data from.
        **kwargs : Any
            |kwargs|
        """
        pass

    def saved_output_files(self, **kwargs: Any) -> list[Path]:
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
        output_file_list = []
        for pattern in self.output_file_pattern:
            output_file_list.extend(self.output_dir.glob(pattern))

        # Add surface image files to the list
        output_file_list.extend(list(self.plot_dir.glob(f"*.{self.output_image_file_extension}")))

        if output_file_list:
            return output_file_list
        else:
            return []

    def output_filename(self, interval: int | None = None, **kwargs: Any) -> str:
        """
        Generate the base output filename for a given interval, but does not include the output directory.

        Parameters
        ----------
        interval : int or None, optional
            The interval number to generate the filename for. If None, the filename will not include an interval number. Default is None.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        str
            The generated output filename for the given interval.
        """
        if self.output_file_prefix is None or self.output_file_extension is None:
            raise ValueError("output_file_prefix and output_file_extension must be set to generate an output filename")
        if interval is not None:
            return f"{self.output_file_prefix}{interval:06d}.{self.output_file_extension}"
        else:
            return f"{self.output_file_prefix}.{self.output_file_extension}"

    def _output_file_processor(self, data_file_list, interval_index: int | None = None, **kwargs: Any) -> xr.Dataset:
        """
        File processor for the output files generated by this component.

        The base method assumes that the output files are in a format that can be read by xarray.open_dataset or xarray.open_mfdataset. Subclasses can override this method to implement custom file processing logic.

        Parameters
        ----------
        interval_index : int or None, optional
            The index of the interval being processed. This can be used to customize the file processor based on the interval. Default is None.
        **kwargs : Any
            |kwargs|
        """
        if interval_index is not None:
            with xr.open_dataset(data_file_list[interval_index]) as ds:
                ds.load()
        elif len(data_file_list) == 1:
            with xr.open_dataset(data_file_list[0]) as ds:
                ds.load()
        else:
            try:
                ds = xr.open_mfdataset(data_file_list, parallel=False, engine="netcdf4", compat="no_conflicts", join="outer")
            except Exception:
                ds = {}  # Fall back when the files cannot be opened together (e.g., because they have different variables or dimensions). In this case, we will open each file separately and concatenate them into a list of Datasets.
                for data_file in tqdm(
                    data_file_list, desc="Reading in files....", unit="files", total=len(data_file_list), position=0, leave=False
                ):
                    with xr.open_mfdataset(data_file, engine="netcdf4") as ds_single:
                        ds[ds_single.interval.item()] = ds_single
        if hasattr(ds, "close"):
            ds.close()
        return ds

    def _overwrite_check(self, files_to_remove: Path | list[Path] | str | list[str]) -> bool:
        """
        Check if the user wants to overwrite existing files or directories, and if they do either from a prompt or based on the `ask_overwrite` attribute, then delete them.

        The user will be prompted to confirm that they want to overwrite existing files or directories if the `ask_overwrite` attribute is set to True. If the user confirms that they want to overwrite the files, then the specified files or directories will be deleted. If the user cancels the operation, then the files will not be deleted and the function will return False. Selecting 'a' will suppress prompts about overwriting files for this operation, but not for future ones.

        Parameters
        ----------
        files_to_remove : Path or list of Path or str or list of str
            The file or list of files that will be removed.

        Returns
        -------
        bool
            True if the user confirmed that they wanted to overwrite the files and that all requested files and/or directories were removed. False if they cancel the operation.
        """
        ask_overwrite = self.ask_overwrite  # This will allow us to temporarily disable prompts about overwriting files for this operation if the user selects 'a' to suppress prompts about overwriting files, without changing the value of self.ask_overwrite for future operations.
        if not isinstance(files_to_remove, list):
            files_to_remove = [files_to_remove]
        for file in files_to_remove:
            file = Path(file)
            if file.exists() and ask_overwrite:
                response = input(
                    f"File '{str(file)}' already exists. To disable this message, set ask_overwrite=False to this instance. Overwrite? (y/N/a): "
                )
                if response.lower() == "a":
                    ask_overwrite = False
                elif response.lower() != "y":
                    print("Operation cancelled by user.")
                    return False
            try:
                if file.is_dir():
                    shutil.rmtree(file)
                else:
                    file.unlink(missing_ok=True)
            except Exception as e:
                print(f"Error removing file {file}: {e}")
                return False
        return True

    def add_save_action(self, action: dict[str, dict]) -> None:
        """
        Add an action to the save_actions property of this component.

        Parameters
        ----------
        action : dict[str, dict]
            A dictionary where the keys are the names of actions that can be performed on this component (e.g. "plot") when calling the save function and the values are the arguments that should be passed to that action when it is called by the save function.
        """
        if not isinstance(action, dict):
            raise TypeError(
                "action must be a dictionary containing a key with a valid action for this component (e.g. 'plot') and the values are the arguments"
            )
        for action_name, args in action.items():
            if hasattr(self, action_name):
                action_method = getattr(self, action_name)
                if not callable(action_method):
                    raise ValueError(f"{action_name} is not a valid action for {self.name}")
            else:
                raise ValueError(f"{action_name} is not a valid action for {self.name}")
            if not isinstance(args, dict):
                raise TypeError(
                    f"Arguments for action {action_name} must be a dictionary of keyword arguments to be passed to the {self.name}.{action_name}() method."
                )
        self._save_actions.append(action)
        return

    @parameter
    def save_actions(self) -> list[dict[str, dict]]:
        """
        A dictionary where the keys are the names of actions that can be performed on this component (e.g. "plot") when calling the save function and the values are the arguments that should be passed to that action when it is called by the save function.
        """
        return self._save_actions

    @save_actions.setter
    def save_actions(self, value) -> None:
        if not isinstance(value, list):
            raise TypeError(
                "save_actions must be a list where each entry is a dictionary containing a key with a valid action for this component (e.g. 'plot') and the values are the arguments"
            )
        for entry in value:
            if not isinstance(entry, dict):
                raise TypeError(
                    "save_actions must be a list where each entry is a dictionary containing a key with a valid action for this component (e.g. 'plot') and the values are the arguments"
                )
            for action, args in entry.items():
                if hasattr(self, action):
                    action_method = getattr(self, action)
                    if not callable(action_method):
                        raise ValueError(f"{action} is not a valid action for {self.name}")
                else:
                    raise ValueError(f"{action} is not a valid action for {self.name}")
                if not isinstance(args, dict):
                    raise TypeError(
                        f"Arguments for action {action} must be a dictionary of keyword arguments to be passed to the {self.name}.{action}() method."
                    )
        self._save_actions = value

    @parameter
    def ask_overwrite(self):
        """
        Flag to indicate whether the user should be prompted to overwrite any old files or not.
        """
        return self._ask_overwrite

    @ask_overwrite.setter
    def ask_overwrite(self, value):
        if value is None:
            value = True
        if not isinstance(value, bool):
            raise TypeError("ask_overwrite must be a bool")
        self._ask_overwrite = value

    def read_saved_output(self, interval: int | None = None, **kwargs) -> tuple[xr.Dataset, ...] | None:
        """
        Read the saved output files for this component and return them as xarray Datasets.

        Parameters
        ----------
        interval : int or None, optional
            The interval number to read. If None, read all intervals. Default is None.
        **kwargs : Any
            |kwargs|

        Returns
        -------
            tuple[xr.Dataset, ...] or None
                A tuple of xarray Datasets containing the data from the saved output files. The order of the Datasets in the tuple corresponds to the order of the file patterns returned by the output_file_pattern property. If no output files are found, returns None.
        """
        import re

        if len(self.output_file_pattern) == 0 or self.output_dir is None:
            return None
        output_ds = []

        for pattern in self._output_file_pattern:
            data_file_list = list(Path(self.output_dir).glob(pattern))
            if len(data_file_list) == 0:
                output_ds.append(None)
                continue
            interval_numbers = []
            if len(data_file_list) > 1:
                for data_file in data_file_list:
                    n = data_file.name.split(self.output_file_prefix)[-1].split(f".{self.output_file_extension}")[0]
                    interval_numbers.append(int(n))
            else:  # handle the special case where the output file pattern does not include an interval number (e.g. "grid.nc")
                n = re.findall(r"\d+", data_file_list[0].name)
                if n:
                    interval_numbers.append(int(n[-1]))

            if len(interval_numbers) > 1:
                tup = sorted(zip(data_file_list, interval_numbers, strict=True))
                data_file_list, interval_numbers = zip(*tup, strict=True)
            data_file_list = list(data_file_list)

            # map the requested interval to the index of interval_numbers
            if interval is not None and len(interval_numbers) > 0:
                if interval < 0:
                    interval_index = interval
                elif interval in interval_numbers:
                    interval_index = interval_numbers.index(interval)
                else:
                    output_ds.append(xr.Dataset())
                    continue
            else:
                interval_index = None

            output_ds.append(self._output_file_processor(data_file_list=data_file_list, interval_index=interval_index, **kwargs))

        if len(output_ds) == 1:
            return output_ds[0]
        else:
            return tuple(output_ds)

    @parameter
    def simdir(self):
        """
        |simdir|

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
    def export_dir(self) -> Path | None:
        """
        The export directory for a component. If None, the component does not have an export directory set.
        """
        if self._export_dir_name is None:
            return None
        export_dir = self.simdir / self._export_dir_name
        if not export_dir.exists():
            try:
                export_dir.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                raise RuntimeError(f"Could not create export directory at {export_dir}") from e
        return export_dir

    @property
    def plot_dir(self) -> Path:
        """
        The directory for plots.
        """
        if self._output_dir_name is None:
            return None
        plot_dir = self.simdir / f"{self._output_dir_name}_images"
        if not plot_dir.exists():
            try:
                plot_dir.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                raise RuntimeError(f"Could not create output directory at {plot_dir}") from e
        return plot_dir

    @property
    def output_image_file_extension(self) -> str:
        """
        The file extension to use for output images.
        """
        return self._output_image_file_extension

    @output_image_file_extension.setter
    def output_image_file_extension(self, value: str):
        """
        Set the file extension to use for output images.

        Parameters
        ----------
        value : str
            The file extension to use for output images (e.g., "png", "jpg", "tif").
        """
        self._output_image_file_extension = value

    @property
    def output_file_pattern(self) -> list[str]:
        """
        Return a list of one or more file patterns that this component will generate when saved.

        This is used when reading save output files or when restting the components. Some components may generate multiple output files with different patterns (e.g. the observed and emplaced crater files generated by the counting component), so this should return a list of all patterns that are generated by this component. If the component does not generate any files, this should return an empty list.

        Returns
        -------
        list[str]
            A list of file patterns (e.g., ["surface*.nc"]) that this component will generate when saved.
            If the component does not generate any files, return an empty list.
        """
        return self._output_file_pattern

    @property
    def output_file_prefix(self) -> str | None:
        """
        Return the file prefix for this component's output files. This is used to identify the files generated by this component when reading saved output or resetting the component. If the component does not generate any files, this should return None.

        Returns
        -------
        str or None
            The file prefix (e.g., "surface") for this component's output files, or None if the component does not generate any files.
        """
        return self._output_file_prefix

    @property
    def output_file_extension(self) -> str | None:
        """
        Return the file extension for this component's output files. This is used to identify the files generated by this component when reading saved output or resetting the component. If the component does not generate any files, this should return None.

        Returns
        -------
        str or None
            The file extension (e.g., "nc") for this component's output files, or None if the component does not generate any files.
        """
        return self._output_file_extension

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
            ask_overwrite=self.ask_overwrite,
        )

    @property
    def driver_to_extension_map(self) -> dict[str, str]:
        """A dictionary mapping valid export/from_file drivers for this component to the corresponding file extensions that can be used when exporting files with that driver."""
        return self._driver_to_extension_map

    @property
    def valid_drivers(self) -> list[str]:
        """A list of valid export/from_file drivers for this component."""
        return list(self._driver_to_extension_map.keys())

    @property
    def valid_extensions(self) -> list[str]:
        """A list of valid export/from_file file extensions for this component."""
        return list(self._driver_to_extension_map.values())


class ComponentBase(CratermakerBase, ABC):
    """
    Base class for components of the Cratermaker project.

    Defines the common parameters and methods for all components in the Cratermaker project, including the maker class that is used to select the correct component from user arguments.


    """

    _registry: dict[str, type[ComponentBase]] = {}

    def __init__(self, **kwargs: Any) -> None:
        """
        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
        **kwargs : Any
            |kwargs|
        """
        object.__setattr__(self, "_save_actions", {})
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
            |kwargs|

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
        |rng|
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        |rng_seed|
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
        |simdir|

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
