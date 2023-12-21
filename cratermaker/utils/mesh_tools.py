from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
import xarray
import jigsawpy
from jigsawpy.savejig import savejig
import sys
import logging
import subprocess
from netCDF4 import Dataset as NetCDFFile
import collections
import importlib.resources

jigsaw_exe = importlib.resources.files('cratermaker').joinpath('bin').joinpath('jigsaw')
mpas_conversion_tool_exe = importlib.resources.files('cratermaker').joinpath('bin').joinpath('MpasMeshConverter.x')


def build_spherical_mesh(cellWidth, lon, lat, earth_radius,
                         out_filename='base_mesh.nc', 
                         dir='./', logger=None):
    """
    Build an MPAS mesh using JIGSAW with the given cell sizes as a function of
    latitude and longitude.

    The result is a mesh file stored in ``out_filename`` as well as several
    intermediate files: ``mesh.log``, ``mesh-HFUN.msh``, ``mesh.jig``,
    ``mesh-MESH.msh``, ``mesh.msh``, and ``mesh_triangles.nc``.

    Parameters
    ----------
    cellWidth : ndarray
        m x n array of cell width in km

    lon : ndarray
        longitude in degrees (length n and between -180 and 180)

    lat : ndarray
        longitude in degrees (length m and between -90 and 90)

    earth_radius : float
        Earth radius in meters

    out_filename : str, optional
        The file name of the resulting MPAS mesh

    dir : str, optional
        A directory in which a temporary directory will be added with files
        produced during mesh conversion and then deleted upon completion.

    logger : logging.Logger, optional
        A logger for the output if not stdout
    """

    with LoggingContext(__name__, logger=logger) as logger:

        da = xarray.DataArray(cellWidth,
                              dims=['lat', 'lon'],
                              coords={'lat': lat, 'lon': lon},
                              name='cellWidth')
        cw_filename = 'cellWidthVsLatLon.nc'
        da.to_netcdf(cw_filename)

        logger.info('Step 1. Generate mesh with JIGSAW')
        jigsaw_driver(cellWidth, lon, lat, on_sphere=True,
                      earth_radius=earth_radius, logger=logger)

        logger.info('Step 2. Convert triangles from jigsaw format to netcdf')
        jigsaw_to_netcdf(msh_filename='mesh-MESH.msh',
                         output_name='mesh_triangles.nc', on_sphere=True,
                         sphere_radius=earth_radius)

        logger.info('Step 3. Convert from triangles to MPAS mesh')
        args = [mpas_conversion_tool_exe,
                'mesh_triangles.nc',
                out_filename]
        check_call(args=args, logger=logger)


def jigsaw_driver(cellWidth, x, y, on_sphere=True, earth_radius=6371.0e3,
                  geom_points=None, geom_edges=None, logger=None):
    """
    A function for building a jigsaw mesh

    Parameters
    ----------
    cellWidth : ndarray
        The size of each cell in the resulting mesh as a function of space

    x, y : ndarray
        The x and y coordinates of each point in the cellWidth array (lon and
        lat for spherical mesh)

    on_sphere : logical, optional
        Whether this mesh is spherical or planar

    earth_radius : float, optional
        Earth radius in meters

    geom_points : ndarray, optional
        list of point coordinates for bounding polygon for planar mesh

    geom_edges : ndarray, optional
        list of edges between points in geom_points that define the bounding polygon

    logger : logging.Logger, optional
        A logger for the output if not stdout
    """
    # Authors
    # -------
    # Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    # setup files for JIGSAW
    opts = jigsawpy.jigsaw_jig_t()
    opts.geom_file = 'mesh.msh'
    opts.jcfg_file = 'mesh.jig'
    opts.mesh_file = 'mesh-MESH.msh'
    opts.hfun_file = 'mesh-HFUN.msh'

    # save HFUN data to file
    hmat = jigsawpy.jigsaw_msh_t()
    if on_sphere:
       hmat.mshID = 'ELLIPSOID-GRID'
       hmat.xgrid = np.radians(x)
       hmat.ygrid = np.radians(y)
    else:
       hmat.mshID = 'EUCLIDEAN-GRID'
       hmat.xgrid = x
       hmat.ygrid = y
    hmat.value = cellWidth
    jigsawpy.savemsh(opts.hfun_file, hmat)

    # define JIGSAW geometry
    geom = jigsawpy.jigsaw_msh_t()
    if on_sphere:
       geom.mshID = 'ELLIPSOID-MESH'
       geom.radii = earth_radius*1e-3*np.ones(3, float)
    else:
       geom.mshID = 'EUCLIDEAN-MESH'
       geom.vert2 = geom_points
       geom.edge2 = geom_edges
    jigsawpy.savemsh(opts.geom_file, geom)

    # build mesh via JIGSAW!
    opts.hfun_scal = 'absolute'
    opts.hfun_hmax = float("inf")
    opts.hfun_hmin = 0.0
    opts.mesh_dims = +2  # 2-dim. simplexes
    opts.optm_qlim = 0.9375
    opts.verbosity = +1

    savejig(opts.jcfg_file, opts)
    check_call([jigsaw_exe, opts.jcfg_file], logger=logger)


def check_call(args, logger=None, log_command=True, timeout=None, **kwargs):
    """
    Call the given subprocess with logging to the given logger.

    Parameters
    ----------
    args : list or str
        A list or string of argument to the subprocess.  If ``args`` is a
        string, you must pass ``shell=True`` as one of the ``kwargs``.

    logger : logging.Logger, optional
        The logger to write output to

    log_command : bool, optional
        Whether to write the command that is running ot the logger

    timeout : int, optional
        A timeout in seconds for the call

    **kwargs : dict
        Keyword arguments to pass to subprocess.Popen

    Raises
    ------
    subprocess.CalledProcessError
        If the given subprocess exists with nonzero status

    """

    if isinstance(args, str):
        print_args = args
    else:
        print_args = ' '.join(args)

    # make a logger if there isn't already one
    with LoggingContext(print_args, logger=logger) as logger:
        if log_command:
            logger.info(f'Running: {print_args}')

        process = subprocess.Popen(args, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, **kwargs)
        stdout, stderr = process.communicate(timeout=timeout)

        if stdout:
            stdout = stdout.decode('utf-8')
            for line in stdout.split('\n'):
                logger.info(line)
        if stderr:
            stderr = stderr.decode('utf-8')
            for line in stderr.split('\n'):
                logger.error(line)

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode,
                                                print_args)


class LoggingContext(object):

    """
    A context manager for creating a logger or using an existing logger

    Attributes
    ----------
    logger : logging.Logger
        A logger that sends output to a log file or stdout/stderr
    """

    def __init__(self, name, logger=None, log_filename=None):
        """
        If ``logger`` is ``None``, create a new logger either to a log file
        or stdout/stderr.  If ``logger`` is anything else, just set the logger
        attribute

        Parameters
        ----------
        name : str
            A unique name for the logger (e.g. ``__name__`` of the calling
            module)

        logger : logging.Logger, optional
           An existing logger that sends output to a log file or stdout/stderr
           to be used in this context

        log_filename : str, optional
            The name of a file where output should be written.  If none is
            supplied, output goes to stdout/stderr
        """
        self.logger = logger
        self.name = name
        self.log_filename = log_filename
        self.handler = None
        self.old_stdout = None
        self.old_stderr = None
        self.existing_logger = logger is not None

    def __enter__(self):
        if not self.existing_logger:
            if self.log_filename is not None:
                # redirect output to a log file
                logger = logging.getLogger(self.name)
                handler = logging.FileHandler(self.log_filename)
            else:
                logger = logging.getLogger(self.name)
                handler = logging.StreamHandler(sys.stdout)

            formatter = MpasFormatter()
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
            logger.propagate = False
            self.logger = logger
            self.handler = handler

            if self.log_filename is not None:
                self.old_stdout = sys.stdout
                self.old_stderr = sys.stderr
                sys.stdout = StreamToLogger(logger, logging.INFO)
                sys.stderr = StreamToLogger(logger, logging.ERROR)
        return self.logger

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.existing_logger:
            if self.old_stdout is not None:
                self.handler.close()
                # restore stdout and stderr
                sys.stdout = self.old_stdout
                sys.stderr = self.old_stderr

            # remove the handlers from the logger (probably only necessary if
            # writeLogFile==False)
            self.logger.handlers = []

        self.stdout = self.original_stdout = sys.stdout
        self.stderr = self.original_stderr = sys.stderr


class MpasFormatter(logging.Formatter):
    """
    A custom formatter for logging
    Modified from:
    https://stackoverflow.com/a/8349076/7728169
    https://stackoverflow.com/a/14859558/7728169
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    # printing error messages without a prefix because they are sometimes
    # errors and sometimes only warnings sent to stderr
    dbg_fmt = "DEBUG: %(module)s: %(lineno)d: %(msg)s"
    info_fmt = "%(msg)s"
    err_fmt = info_fmt

    def __init__(self, fmt=info_fmt):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._fmt = MpasFormatter.dbg_fmt

        elif record.levelno == logging.INFO:
            self._fmt = MpasFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._fmt = MpasFormatter.err_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result


class StreamToLogger(object):
    """
    Modified based on code by:
    https://www.electricmonk.nl/log/2011/08/14/redirect-stdout-and-stderr-to-a-logger-in-python/
    Copyright (C) 2011 Ferry Boender
    License: GPL, see https://www.electricmonk.nl/log/posting-license/
    Fake file-like stream object that redirects writes to a logger instance.
    """

    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, str(line.rstrip()))

    def flush(self):
        pass



def jigsaw_to_netcdf(msh_filename, output_name, on_sphere, sphere_radius=None):
    """
    Converts mesh data defined in triangle format to NetCDF

    Parameters
    ----------
    msh_filename : str
        A JIGSAW mesh file name
    output_name: str
        The name of the output file
    on_sphere : bool
        Whether the mesh is spherical or planar
    sphere_radius : float, optional
        The radius of the sphere in meters.  If ``on_sphere=True`` this argument
        is required, otherwise it is ignored.
    """
    # Authors: Phillip J. Wolfram, Matthew Hoffman and Xylar Asay-Davis

    grid = NetCDFFile(output_name, 'w', format='NETCDF3_CLASSIC')

    # Get dimensions
    # Get nCells
    msh = readmsh(msh_filename)
    nCells = msh['POINT'].shape[0]

    # Get vertexDegree and nVertices
    vertexDegree = 3  # always triangles with JIGSAW output
    nVertices = msh['TRIA3'].shape[0]

    if vertexDegree != 3:
        ValueError("This script can only compute vertices with triangular "
                   "dual meshes currently.")

    grid.createDimension('nCells', nCells)
    grid.createDimension('nVertices', nVertices)
    grid.createDimension('vertexDegree', vertexDegree)

    # Create cell variables and sphere_radius
    xCell_full = msh['POINT'][:, 0]
    yCell_full = msh['POINT'][:, 1]
    zCell_full = msh['POINT'][:, 2]
    for cells in [xCell_full, yCell_full, zCell_full]:
        assert cells.shape[0] == nCells, 'Number of anticipated nodes is' \
                                         ' not correct!'
    if on_sphere:
        grid.on_a_sphere = "YES"
        grid.sphere_radius = sphere_radius
        # convert from km to meters
        xCell_full *= 1e3
        yCell_full *= 1e3
        zCell_full *= 1e3
    else:
        grid.on_a_sphere = "NO"
        grid.sphere_radius = 0.0

    # Create cellsOnVertex
    cellsOnVertex_full = msh['TRIA3'][:, :3] + 1
    assert cellsOnVertex_full.shape == (nVertices, vertexDegree), \
        'cellsOnVertex_full is not the right shape!'

    # Create vertex variables
    xVertex_full = np.zeros((nVertices,))
    yVertex_full = np.zeros((nVertices,))
    zVertex_full = np.zeros((nVertices,))

    for iVertex in np.arange(0, nVertices):
        cell1 = cellsOnVertex_full[iVertex, 0]
        cell2 = cellsOnVertex_full[iVertex, 1]
        cell3 = cellsOnVertex_full[iVertex, 2]

        x1 = xCell_full[cell1 - 1]
        y1 = yCell_full[cell1 - 1]
        z1 = zCell_full[cell1 - 1]
        x2 = xCell_full[cell2 - 1]
        y2 = yCell_full[cell2 - 1]
        z2 = zCell_full[cell2 - 1]
        x3 = xCell_full[cell3 - 1]
        y3 = yCell_full[cell3 - 1]
        z3 = zCell_full[cell3 - 1]

        pv = circumcenter(on_sphere, x1, y1, z1, x2, y2, z2, x3, y3, z3)
        xVertex_full[iVertex] = pv.x
        yVertex_full[iVertex] = pv.y
        zVertex_full[iVertex] = pv.z

    meshDensity_full = grid.createVariable(
        'meshDensity', 'f8', ('nCells',))

    for iCell in np.arange(0, nCells):
        meshDensity_full[iCell] = 1.0

    del meshDensity_full

    var = grid.createVariable('xCell', 'f8', ('nCells',))
    var[:] = xCell_full
    var = grid.createVariable('yCell', 'f8', ('nCells',))
    var[:] = yCell_full
    var = grid.createVariable('zCell', 'f8', ('nCells',))
    var[:] = zCell_full
    var = grid.createVariable('xVertex', 'f8', ('nVertices',))
    var[:] = xVertex_full
    var = grid.createVariable('yVertex', 'f8', ('nVertices',))
    var[:] = yVertex_full
    var = grid.createVariable('zVertex', 'f8', ('nVertices',))
    var[:] = zVertex_full
    var = grid.createVariable(
        'cellsOnVertex', 'i4', ('nVertices', 'vertexDegree',))
    var[:] = cellsOnVertex_full

    grid.sync()
    grid.close()
    
    
def readmsh(fname):
    """
    Reads JIGSAW msh structure and produces a dictionary with values.

    Phillip J. Wolfram
    09/22/2017
    """

    dataset = {}
    datavals = {}
    datavals['HEADER'] = ';'
    datavals['ARRAY'] = None
    with open(fname) as f:
        line = f.readline()
        while line:
            if line[0] == '#':
                datavals['HEADER'] += line[1:] + ';'
                line = f.readline()
                continue
            if '=' in line:
                datavals, dataset = _store_datavals(datavals, dataset)
                if 'COORD' in line:
                    name = 'COORD' + line.split('=')[1][0]
                    datavals[name] = line.split(';')[-1]
                else:
                    vals = line.split('=')
                    value = vals[1] if ';' in vals[1] else int(vals[1])
                    datavals[vals[0]] = value
                line = f.readline()
                continue

            # just numbers
            arrayvals = np.asarray(line.split(';'), dtype='f8')
            if datavals['ARRAY'] is None:
                datavals['ARRAY'] = [arrayvals]
            else:
                datavals['ARRAY'].append(arrayvals)
            line = f.readline()
            continue
        datavals, dataset = _store_datavals(datavals, dataset)

    return dataset


def _store_datavals(datavals, dataset):  # {{{

    if datavals['ARRAY'] is not None:
        # remove empty data
        if np.all(datavals['ARRAY'] == np.array(None, dtype='object')):
            datavals.pop('ARRAY')
        for key in [aval for aval in datavals.keys()
                    if aval in ['HEADER', 'MSHID', 'NDIMS']]:
            if key in dataset:
                dataset[key] += datavals[key]
            else:
                dataset[key] = datavals[key]
            datavals.pop(key)
        entryname = [aval for aval in datavals.keys() if aval not in [
            'ARRAY']]

        if 'TRI' in entryname[0]:
            dtype = 'i'
        else:
            dtype = 'f8'
        datavals['ARRAY'] = np.asarray(datavals['ARRAY'], dtype=dtype)

        # decided to throw away "index" from msh because it isn't truly a
        # real number
        dataset[entryname[0]] = datavals['ARRAY']
        datavals = {}
        datavals['ARRAY'] = None

    return datavals, dataset  # }}}



point = collections.namedtuple('Point', ['x', 'y', 'z'])


def circumcenter(on_sphere, x1, y1, z1, x2, y2, z2, x3, y3, z3):
    """
    Compute the circumcenter of the triangle (possibly on a sphere)
    with the three given vertices in Cartesian coordinates.

    Returns
    -------
    center : point
        The circumcenter of the triangle with x, y and z attributes

    """
    p1 = point(x1, y1, z1)
    p2 = point(x2, y2, z2)
    p3 = point(x3, y3, z3)
    if on_sphere:
        a = (p2.x - p3.x)**2 + (p2.y - p3.y)**2 + (p2.z - p3.z)**2
        b = (p3.x - p1.x)**2 + (p3.y - p1.y)**2 + (p3.z - p1.z)**2
        c = (p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2

        pbc = a * (-a + b + c)
        apc = b * (a - b + c)
        abp = c * (a + b - c)

        xv = (pbc * p1.x + apc * p2.x + abp * p3.x) / (pbc + apc + abp)
        yv = (pbc * p1.y + apc * p2.y + abp * p3.y) / (pbc + apc + abp)
        zv = (pbc * p1.z + apc * p2.z + abp * p3.z) / (pbc + apc + abp)
    else:
        d = 2 * (p1.x * (p2.y - p3.y) + p2.x *
                 (p3.y - p1.y) + p3.x * (p1.y - p2.y))

        xv = ((p1.x**2 + p1.y**2) * (p2.y - p3.y) + (p2.x**2 + p2.y**2)
              * (p3.y - p1.y) + (p3.x**2 + p3.y**2) * (p1.y - p2.y)) / d
        yv = ((p1.x**2 + p1.y**2) * (p3.x - p2.x) + (p2.x**2 + p2.y**2)
              * (p1.x - p3.x) + (p3.x**2 + p3.y**2) * (p2.x - p1.x)) / d
        zv = 0.0

        # Optional method to use barycenter instead.
        # xv = p1.x + p2.x + p3.x
        # xv = xv / 3.0
        # yv = p1.y + p2.y + p3.y
        # yv = yv / 3.0
    return point(xv, yv, zv)

