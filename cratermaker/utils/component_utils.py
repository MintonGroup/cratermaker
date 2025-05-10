from __future__ import annotations

import importlib
import pkgutil
from abc import ABC
from typing import Any

from cratermaker.core.base import CratermakerBase
from cratermaker.utils.general_utils import parameter


class ComponentBase(CratermakerBase, ABC):
    _registry: dict[str, type[ComponentBase]] = {}

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)

    def __str__(self) -> str:
        # Get the name of the class
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
                raise KeyError(
                    f"Unknown component model: {component}. Available models: {cls.available()}"
                )
            return cls._registry[component](**kwargs)
        elif isinstance(component, type) and issubclass(component, ComponentBase):
            return component(**kwargs)
        elif isinstance(component, ComponentBase):
            return component
        else:
            raise TypeError(
                f"component must be a string or a subclass of component, not {type(component)}"
            )

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


def import_components(
    package_name: str, package_path: list[str], ignore_private: bool = True
) -> None:
    """
    Import all modules of a component package, optionally skipping private modules.

    Parameters
    ----------
    package_name : str
        The full name of the package (e.g., "cratermaker.components.component").
    package_path : list[str]
        The __path__ attribute of the package (usually just __path__).
    ignore_private : bool, optional
        Whether to ignore modules whose names start with an underscore (default is True).
    """
    for _, module_name, _ in pkgutil.iter_modules(package_path):
        if ignore_private and module_name.startswith("_"):
            continue
        importlib.import_module(f"{package_name}.{module_name}")
