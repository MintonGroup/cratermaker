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

# Map of OGR drivers to file extensions
VECTOR_DRIVER_TO_EXTENSION_MAP = {
    "PCIDSK": "pix",
    "PDS4": "xml",
    "PDF": "pdf",
    "MBTILES": "mbtiles",
    "ESRI SHAPEFILE": "shp",
    "MAPINFO FILE": "tab",
    "S57": "000",
    "DGN": "dgn",
    "CSV": "csv",
    "GML": "gml",
    "GPX": "gpx",
    "KML": "kml",
    "GEOJSON": "json",
    "GEOJSONSEQ": "geojsonl",
    "OGR_GMT": "gmt",
    "GPKG": "gpkg",
    "SQLite": "sqlite",
    "WASP": "map",
    "OPENFILEGDB": "gdb",
    "DXF": "dxf",
    "FLATGEOBUF": "fgb",
    "PGDUMP": "sql",
    "GPSBABEL": "mps",
    "ODS": "ods",
    "XLSX": "xlsx",
    "JML": "jml",
    "VDV": "txt",
    "MVT": "mvt",
    "PMTiles": "pmtiles",
    "JSONFG": "json",
    "MIRAMONVECTOR": "pol",
}

PYVISTA_SHOW_KWARGS = [
    "title",
    "window_size",
    "interactive",
    "auto_close",
    "interactive_update",
    "full_screen",
    "screenshot",
    "return_img",
    "cpos",
    "jupyter_backend",
    "return_viewer",
    "return_cpos",
    "before_close_callback",
    "store_image_depth",
]

PYVISTA_ADD_MESH_KWARGS = [
    "color",
    "style",
    "scalars",
    "clim",
    "show_edges",
    "edge_color",
    "point_size",
    "line_width",
    "opacity",
    "flip_scalars",
    "lighting",
    "n_colors",
    "interpolate_before_map",
    "cmap",
    "label",
    "reset_camera",
    "scalar_bar_args",
    "show_scalar_bar",
    "multi_colors",
    "name",
    "texture",
    "render_points_as_spheres",
    "render_lines_as_tubes",
    "smooth_shading",
    "split_sharp_edges",
    "ambient",
    "diffuse",
    "specular",
    "specular_power",
    "nan_color",
    "nan_opacity",
    "culling",
    "rgb",
    "categories",
    "silhouette",
    "use_transparency",
    "below_color",
    "above_color",
    "annotationspickable",
    "preference",
    "log_scale",
    "pbr",
    "metallic",
    "roughness",
    "render",
    "user_matrix",
    "component",
    "emissive",
    "copy_mesh",
    "backface_params",
    "show_vertices",
    "edge_opacity",
    "remove_existing_actor",
]


# Optional: controlled public API
__all__ = []
