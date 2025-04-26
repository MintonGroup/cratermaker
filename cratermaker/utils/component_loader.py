import pkgutil
import importlib

def import_components(package_name: str, package_path: list[str], ignore_private: bool = True) -> None:
    """
    Import all modules of a component package, optionally skipping private modules.

    Parameters
    ----------
    package_name : str
        The full name of the package (e.g., "cratermaker.components.morphology").
    package_path : list[str]
        The __path__ attribute of the package (usually just __path__).
    ignore_private : bool, optional
        Whether to ignore modules whose names start with an underscore (default is True).
    """
    for _, module_name, _ in pkgutil.iter_modules(package_path):
        if ignore_private and module_name.startswith("_"): 
            continue
        importlib.import_module(f"{package_name}.{module_name}")