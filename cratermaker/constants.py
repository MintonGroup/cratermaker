# Default filenames and paths 
_CONFIG_FILE_NAME = "cratermaker.yaml"
_DATA_DIR = "surface_data"
_EXPORT_DIR = "export"
_COMBINED_DATA_FILE_NAME = "surface.nc"
_GRID_FILE_NAME = "grid.nc"
_CIRCLE_FILE_NAME = "circles.vtp"
_VTK_FILE_EXTENSION = "vtp"
_COMPONENT_NAMES = ["target", "scaling", "production", "morphology", "projectile", "surface"]

_SMALLFAC = 1.0e-5 # This is a factor used to determine the smallest length scale in the grid

# Optional: controlled public API
__all__ = []