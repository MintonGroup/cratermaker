import numpy as np

# Custom type aliases for better readability
FloatLike = float | int | np.number
PairOfFloats = tuple[float, float] | list[float] | np.ndarray

# Default filenames and paths
_CONFIG_FILE_NAME = "cratermaker.yaml"
_COMPONENT_NAMES = [
    "counting",
    "target",
    "scaling",
    "production",
    "morphology",
    "projectile",
    "surface",
]
# This is a factor used to determine the smallest length scale in the grid
_SMALLFAC = 1.0e-5
_VSMALL = 10 * np.finfo(np.float64).tiny
_LOGVSMALL = np.log10(_VSMALL)


# Optional: controlled public API
__all__ = []
