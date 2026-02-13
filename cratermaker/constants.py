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
EXPORT_DRIVER_TO_EXTENSION_MAP = {
    "PCIDSK": "pix",
    "PDS4": "xml",
    "PDF": "pdf",
    "MBTiles": "mbtiles",
    "ESRI Shapefile": "shp",
    "MapInfo File": "tab",
    "S57": "000",
    "DGN": "dgn",
    "CSV": "csv",
    "GML": "gml",
    "GPX": "gpx",
    "KML": "kml",
    "GeoJSON": "json",
    "GeoJSONSeq": "geojsonl",
    "OGR_GMT": "gmt",
    "GPKG": "gpkg",
    "SQLite": "sqlite",
    "WAsP": "map",
    "OpenFileGDB": "gdb",
    "DXF": "dxf",
    "FlatGeobuf": "fgb",
    "PGDUMP": "sql",
    "GPSBabel": "mps",
    "ODS": "ods",
    "XLSX": "xlsx",
    "JML": "jml",
    "VDV": "txt",
    "MVT": "mvt",
    "PMTiles": "pmtiles",
    "JSONFG": "json",
    "MiraMonVector": "pol",
}


# Optional: controlled public API
__all__ = []
