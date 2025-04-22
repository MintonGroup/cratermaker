import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.spatial.transform import Rotation as R
from typing import Any
from numpy.typing import NDArray
from cratermaker.utils.general_utils import validate_and_convert_location, parameter
from cratermaker.utils.custom_types import FloatLike, PairOfFloats
from cratermaker.components.grid import register_grid_type, GridMaker

@register_grid_type("hireslocal")
class HiResLocalGrid(GridMaker):
    """
    Create a uniform grid configuration with the given pixel size.
    
    Parameters
    ----------
    pix : FloatLike
        The approximate face size inside the local region in meters.
    radius: FloatLike
        The radius of the target body in meters.
    local_radius : FloatLike
        The radius of the local region in meters.
    local_location : PairOfFloats
        The longitude and latitude of the location in degrees.
    superdomain_scale_factor : FloatLike
        A factor defining the ratio of cell size to the distance from the local boundary. This is set so that smallest craters 
        that are modeled outside the local region are those whose ejecta could just reach the boundary.
        
    Returns
    -------
    HiResLocalGrid
        An instance of the HiResLocalGrid clas initialized with the given point distribution
    """        
    def __init__(self, 
                 pix: FloatLike, 
                 radius: FloatLike, 
                 local_radius: FloatLike, 
                 local_location: PairOfFloats,
                 superdomain_scale_factor: FloatLike,
                 **kwargs: Any):
        super().__init__(**kwargs)
        self.pix = pix
        self.radius = radius        
        self.local_radius = local_radius
        self.local_location = local_location
        self.superdomain_scale_factor = superdomain_scale_factor
        
    def _generate_variable_size_array(self) -> tuple[NDArray, NDArray, NDArray]:
        """
        Create an array of target pixel sizes for pairs of longitude and latitude values for a high-resolution local mesh around a 
        given location, with adaptive cell sizes outside the local region.
        
            
        Returns
        -------
        pix_array: ndarray
            m x n array of pixel sizes
        lon : ndarray
            longitude in degrees (length n and between -180 and 180)
        lat : ndarray
            latitude in degrees (length m and between -90 and 90)
        """
        
        def _pix_func(lon,lat):
            from cratermaker.core.surface import Surface
            lon_rad = np.radians(lon)
            lat_rad = np.radians(lat)
            loc_lon_rad = 0.0 #np.radians(self.local_location[0])
            loc_lat_rad = 0.0 #np.radians(self.local_location[1])

            # Calculate distance from the location to the grid point
            distance = Surface.calculate_haversine_distance(loc_lon_rad, loc_lat_rad, lon_rad, lat_rad, self.radius)
            ans = np.where(distance <= self.local_radius, self.pix, (distance - self.local_radius) / self.superdomain_scale_factor + self.pix)
            return ans
        

        # Suppose we know pix(lat, lon)
        # Step 1: Construct a fine preliminary grid to estimate integrals
       
        Lat = np.linspace(-90., 90., 733)
        Lon = np.linspace(-180., 180., 733)
        LAT, LON = np.meshgrid(Lat, Lon, indexing='ij')

        pix_values = _pix_func(LON, LAT)   # Evaluate pix on this fine grid
        w = 1.0 / pix_values

        # Step 2: Integrate w over longitude to get W_lat(lat)
        W_lat_vals = np.trapezoid(w, x=Lon, axis=1)  # integrate along lon dimension
        W_lat_cumulative = np.cumsum(W_lat_vals)
        W_lat_cumulative -= W_lat_cumulative[0] # normalize from 0 to 1
        W_lat_cumulative /= W_lat_cumulative[-1]  

        # Create a function to invert W_lat(lat)
        f_lat = interp1d(W_lat_cumulative, Lat, bounds_error=False, fill_value='extrapolate')
        
        M = int(2*np.pi*self.radius/pix_values.max()) - 1
        while M > 0:
            badval=False

            # Step 3: For each lat interval, choose lon lines similarly
            N = M + 1  # number of lon lines
            lat_lines = f_lat(np.linspace(0, 1, N))
            lon_lines = np.zeros((M, N))

            for i in range(M):
                lat_low, lat_high = lat_lines[i], lat_lines[i+1]
                # Extract w in this lat band
                mask = (LAT >= lat_low) & (LAT <= lat_high)
                w_band = w[mask].reshape(-1, len(Lon))  
                # Integrate this band over lat
                w_band_vals = np.trapezoid(w_band, x=Lat[(Lat>=lat_low)&(Lat<=lat_high)], axis=0)
                W_lon_band_cumulative = np.cumsum(w_band_vals)
                if W_lon_band_cumulative[-1] > 0:
                    W_lon_band_cumulative -= W_lon_band_cumulative[0]
                    W_lon_band_cumulative /= W_lon_band_cumulative[-1]
                    f_lon = interp1d(W_lon_band_cumulative, Lon, bounds_error=False, fill_value='extrapolate')
                    lon_lines[i,:] = f_lon(np.linspace(0, 1, N))
                else:
                    badval=True
                    M = M // 2
                    break
                
            if not badval:
                break
            
        if badval:
            raise ValueError("Could not generate a grid with the given parameters. Please try again with different parameters.")
            
        LAT = np.zeros((M, N))
        LON = np.zeros((M, N))
        
        for i in range(M):
            LAT[i, :] = lat_lines[i]  # every point on this horizontal line has the same lat

        LON = lon_lines
            
        pix_array = _pix_func(LON, LAT) 

        return pix_array, LON, LAT

    def _rotate_point_cloud(self,points):
        """
        Rotate a point cloud so that the point at [-r,0,0] moves to (lon,lat)
        using the convention:
        - Longitude [-180, 180] degrees, increasing eastward
        - Latitude [-90, 90] degrees, increasing northward
        - (0,0) lon/lat corresponds to [r,0,0] in Cartesian space
        
        Parameters:
            points (np.ndarray): Nx3 array of (x,y,z) points.
            lon (float): Target longitude in degrees.
            lat (float): Target latitude in degrees.
            r (float): Radius of the sphere (default=1).
        
        Returns:
            np.ndarray: Rotated Nx3 point cloud.
        """

        # Convert target lon, lat to radians
        lon_rad, lat_rad = np.radians(self.local_location)
        
        # Compute target unit vector (correcting for lon,lat convention)
        target = np.array([
            self.radius * np.cos(lat_rad) * np.cos(lon_rad),  
            self.radius * np.cos(lat_rad) * np.sin(lon_rad),  
            self.radius * np.sin(lat_rad)                     
        ])

        # Original vector (the point we want to move)
        original = np.array([-self.radius, 0, 0])  # Starts at [-r, 0, 0]
        if np.isclose(lon_rad, 0.0) and np.isclose(lat_rad, 0.0):
            rotation = R.from_euler('z', 180, degrees=True)  # 180-degree rotation around z-axis
            return rotation.apply(points)

        # Compute the axis of rotation (cross product)
        axis = np.cross(original, target)
        axis_norm = np.linalg.norm(axis)
        
        # If the axis is zero (no rotation needed), return the original points
        if axis_norm < 1e-10:
            return points

        axis /= axis_norm  # Normalize the axis

        # Compute the rotation angle (dot product)
        angle = np.arccos(np.clip(np.dot(original / self.radius, target / self.radius), -1.0, 1.0))  # Normalize for dot product

        # Create rotation object
        rotvec = axis * angle  # Convert to rotation vector
        rotation = R.from_rotvec(rotvec)  # Create rotation from axis-angle

        return rotation.apply(points)

    def generate_face_distribution(self, **kwargs: Any) -> NDArray:
        """
        Creates the points that define the mesh centers.
           
        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.
        
        """                
        
        print(f"Generating a mesh with variable resolution faces")
        print(f"Center of local region: {self.local_location}")
        print(f"Size of local region: {self.local_radius:.2f} m")
        print(f"Hires region pixel size: {self.pix:.2f} m")
        print(f"Lores region pixel size: {self.pix * self.superdomain_scale_factor:.2f} m")
        pix_array, lon, lat = self._generate_variable_size_array()
        
        points = []
        n = lon.shape[0]
        m = lat.shape[1]
       
        for i in range(n-1):
            for j in range(m-1):
                lon_range = (lon[i,j], lon[i,j+1])
                lat_range = (lat[i,j], lat[i+1,j])
                p = self._distribute_points(distance=pix_array[i,j]/self.radius, lon_range=lon_range, lat_range=lat_range) 
                if p is not None:
                    points.append(p)
                
        points = np.concatenate(points, axis=1)
        points = self._rotate_point_cloud(points.T).T
        
        return points

    def generate_grid(self,
                      grid_file: os.PathLike,
                      grid_hash: str | None = None,
                      **kwargs: Any) -> tuple[os.PathLike, os.PathLike]:        
        super().generate_grid(grid_file=grid_file, grid_hash=grid_hash, **kwargs)
        face_areas = self.grid.face_areas 
        face_sizes = np.sqrt(face_areas / (4 * np.pi))
        pix_min = face_sizes.min().item() * self.radius
        pix_max = face_sizes.max().item() * self.radius
        print(f"Generated {self.grid.n_face} faces")
        print(f"Effective pixel size range: {pix_min:.2f},{pix_max:.2f} m")
        return

    @parameter
    def pix(self):
        """
        The approximate face size for a cell of the mesh.
        """
        return self._pix
    
    @pix.setter
    def pix(self, value: FloatLike):
        if not isinstance(value, FloatLike) or np.isnan(value) or np.isinf(value) or value <= 0:
            raise TypeError("pix must be a positive float")
        self._pix = value

    @parameter
    def local_radius(self):
        """
        The radius of the local region in meters.
        """
        return self._local_radius

    @local_radius.setter
    def local_radius(self, value: FloatLike):
        if not isinstance(value, FloatLike) or np.isnan(value) or np.isinf(value) or value <= 0:
            raise TypeError("local_radius must be a positive float")
        if value > 2 * np.pi * self.radius:
            raise ValueError("local_radius must be less than 2 * pi * radius of the target body")
        self._local_radius = value
        
    @parameter
    def local_location(self):
        """
        The longitude and latitude of the location in degrees.
        """
        return self._local_location
    
    @local_location.setter
    def local_location(self, value: PairOfFloats):
        if not isinstance(value, tuple) or len(value) != 2:
            raise TypeError("local_location must be a tuple of two floats")
        self._local_location = validate_and_convert_location(value)
        
    @parameter
    def superdomain_scale_factor(self):
        """
        A factor defining the ratio of cell size to the distance from the local boundary. This is set so that smallest craters that are 
        modeled outside the local region are those whose ejecta could just reach the boundary.
        """
        return self._superdomain_scale_factor
    
    @superdomain_scale_factor.setter
    def superdomain_scale_factor(self, value: FloatLike):
        if not isinstance(value, FloatLike) or np.isnan(value) or np.isinf(value) or value < 1.0:
            raise TypeError("superdomain_scale_factor must be a positive float greater than or equal to 1")
        self._superdomain_scale_factor = value
   