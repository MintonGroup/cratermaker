.. currentmodule:: cratermaker

.. image:: ../_static/surface_grid.png
    :alt: Surface
    :align: center
    :width: 300px

.. _ug-surface:

Surface
=======

Cratermaker's :ref:`Surface <api-surface>` component is used to represent target body's topography and other properties of its surface. Its prupose is to handle surface-related data by providing methods for setting elevation data, and calculating surface-related questions. This tool contains three classes of surface implementations that can be chosen by the user: **icosphere**, **arbitrary_resolution**, **hireslocal**. 

UxArray-based surface mesh
--------------------------

The surface of a celestial body in Cratermaker is represented as a sphere that has been discretized as an unstructured polygonal mesh using the `UxArray <https://uxarray.readthedocs.io/en/latest/index.html>`_ package. UxArray provides a rich set of tools for representing unstructured mesh geometry and data associated with the mesh, through their `UxDataset <https://uxarray.readthedocs.io/en/latest/user-guide/data-structures.html#uxdataset>`_ and associated `Grid <https://uxarray.readthedocs.io/en/latest/user-guide/data-structures.html#grid>`_.  The surface mesh is composed of faces, nodes, and edges, where each face is a polygonal shape defined by its nodes. The faces are connected to each other through edges, and the nodes are the points in space that define the corners of the faces. A simple diagram showing the relationship between faces, nodes, and edges is shown below:

.. image:: ../_static/mesh_diagram.svg
    :alt: Surface faces, nodes, and edges
    :align: center
    :width: 300px

In the above image, show a single face with 6 nodes and 6 edges, surrounded by 6 neighboring faces. Each face, node, and edge is identified with an integer index. A Surface object contains a number of attributes that represent the mesh geometry and connectivity, such as:

- :attr:`n_face`: The number of faces in the surface mesh.
- :attr:`n_node`: The number of nodes in the surface mesh.
- :attr:`n_edge`: The number of edges in the surface mesh.
- :attr:`face_area` : An array that contains the area of each face in m\ :sup:`2`
- :attr:`face_size` : The "effective" size of each face, which is defined as the square root of the area of each face in meters.
- :attr:`face_x`, :attr:`face_y`, :attr:`face_z`: Arrays that contain the x, y, and z coordinates of the center of each face in Cartesian coordinates.
- :attr:`node_x`, :attr:`node_y`, :attr:`node_z`: Arrays that contain the x, y, and z coordinates of each node in Cartesian coordinates.
- :attr:`face_lon`, :attr:`face_lat`: Arrays that contain the longitude and latitude of the center of each face in degrees.
- :attr:`node_lon`, :attr:`node_lat`: Arrays that contain the longitude and latitude of each node in degrees.
- :attr:`face_elevation`, :attr:`node_elevation`: Arrays that contain the elevation of each face and node in meters. The mesh itself remains static inside Cratermaker, and so the elevations are only applied to the mesh when it is visualized.
- :attr:`edge_length`: An array that contains the length of each edge in meters.
- :attr:`edge_face_distance`: An array that contains the distance between the centers of the two faces that saddle each edge in meters.
- :attr:`face_face_connectivity`: An array that contains the faces that surround each face.
- :attr:`face_node_connectivity`: An array that contains the nodes that are associated with each face.  
- :attr:`node_face_connectivity`: An array that contains the faces that are associated with each node.

The default Surface type: Icosphere
-----------------------------------
The default Surface is "icosphere" with a grid level of 8 and the :ref:`Target <ug-target>` being the Moon. The default can be set by: 

.. ipython:: python
    :okwarning:

    from cratermaker import Surface
    surface=Surface.maker()
    print(surface)

This is equivalent to:

.. ipython:: python
    :okwarning:

    surface=Surface.maker(surface='icosphere', gridlevel=8, target='Moon')
    print(surface)

This surface consists of a uniform grid configuration with polygonal faces that will be subdivided by the input value for the gridlevel argument.  Though it is limited to a few resolutions, the icosphere surface will have the most uniform face sizes. The number of faces and nodes of the icosphere is determined by the formulas :math:`N_{face} = 10\times4^{gridlevel}+2` and :math:`N_{node} = 20\times4^{gridlevel}`. The Surface object contains an attribute called `pix`, which is the value of the "effective pixel size" in meters, where :math:`pix=\sqrt{\left<Area_{face}\right>}`. The following table shows the relationship between the grid level, the number of faces and nodes, and the effective pixel size for a default target of the Moon

.. csv-table::
    :header: "gridlevel", "n_face", "n_node", "pix (for Moon target)"
    :widths: 10, 10, 10, 10

    5, 10242, 20480, 60.8 km ± 2.42 km
    6, 40962, 81920, 30.4 km ± 1.25 km
    7, 163842, 327680, 15.2 km ± 634 m
    8, 655362, 1310720, 7.60 km ± 319 m
    9, 2621442, 5242880, 3.80 km ± 160 m


Lower values of `gridlevel` will result in fewer but larger face sizes, which can be computed quickly but will not resolve detail well. Higher values of `gridlevel` will result in more faces with smaller areas, which will resolve detail better but will take longer to generate and use, and will consume more memory. We recommend keeping `gridlevel` to between ~7-9. Also, keep in mind that the value of `pix` in the table above is computed for the Moon, and will vary for other targets.


Arbitrary resolution grids
--------------------------

While the `icosphere` surface generates the most regular grids, it is limited to only a few fixed face sizes. If you wish to have more control over the sizes of your faces, you can use the "arbitrary_resolution" surface type instead of "icosphere." The "arbitrary_resolution" surface takes an argument called `pix`, which sets the approximate size of the faces of the grid. The value of `pix` is given in units of meter, and it is defined such that the area of each face will on average be :math:`pix^2`.  For instance, suppose we want to create a surface representation of planet Mercury with a resolution of ~20 km per face:

.. ipython:: python
    :okwarning:

    surface=Surface.maker(surface='arbitrary_resolution', target='Mercury', pix=20e3)
    print(surface)


As seen above, the arbitrary resolution grid also creates a uniform grid configurations but allows the user to define their pixel size. Unlike the icosphere, the faces on the surface will be more irregular in shape, making it less ideal. 

High resolution local grids
-----------------------------

In many application of Cratermaker, it is useful to model a small local region at high resolution. This surface type requires the following 4 arguments:

- `pix`: The pixel size in meters within the high resolution local region.
- `local_radius`: The radius of the local region in meters.
- `local_location`: The center of the local region in degrees, given as a tuple of (longitude, latitude).
- `superdomain_scale_factor`: A scale factor that determines how much larger the faces will be at the antipode of the local region.

.. note::
    When used as part of a :ref:`Simulation <ug-simulation>`, the `superdomain_scale_factor` can be omitted. It is computed using the Simulation's production and morphology models in order to compute the sizes of faces in the superdomain as a function of distance from the local region boundary. 

For instance, suppose we want to generate a high resolution local grid on the Moon with a resolution of ~10 m per pixel, with a radius of 5 km, centered at the equator and prime meridian (0° longitude, 0° latitude), and a superdomain scale factor of 10,000:

.. ipython:: python
    :okwarning:

    surface=Surface.maker("hireslocal", pix=10.0, local_radius=5e3, local_location=(0, 0), superdomain_scale_factor=10000)
    print(surface)

The "hireslocal" surface type works somewhat differently than the others. For instance, the diffusive degradation is only applied on the local region. You can think of the local region as the "primary" surface being modeled, and the superdomain as simply a source for distal ejecta fram large far away craters. 


Extracting a local subsection of the surface
--------------------------------------------

Many of the operations that Cratermaker does on a surface only affect a small portion of the full grid at a time. The Surface class has a tool that is used to efficiently extract a local subsection of a given surface without making a copy in memory. This is done by creating a :class:`LocalSurface` object, which is a view of the original surface that only contains the faces and nodes within a specified radius of a given location. To generate a LocalSurface, you can use the :meth:`extract_region` method of the Surface class. This method takes two arguments: `location`, which is a tuple of (longitude, latitude) in degrees, and `region_radius`, which is the radius of the region in meters.

For instance, suppose we'd like to extract a 16 km radius region at the south pole of the Moon:

.. ipython:: python
    :okwarning:

    surface=Surface.maker()
    print(surface)
    region=surface.extract_region(location=(0,-90),region_radius=16e3)
    print(region)

As we can see, this selects only 15 of the full 655362 faces, which is a significant reduction in the number of faces and nodes that need to be processed. The LocalSurface object can be used to perform operations on this local region without affecting the full surface. 

.. note::
    The "hireslocal" Surface type contains a built-in attribute called `local`, which represents the high resolution region of the surface. In addition, when `extract_region` is called on a "hireslocal" surface, it will return a special `LocalHiResLocalSurface` object that contains within it an additional object called `local_overlap`. This is a view of only the portion of the extracted region that overlaps the high resolution region (or None if there is no overlap).

Using a Surface object
----------------------
Once you have either a Surface or LocalSurface object, you are now able to perform numerous surface-related computations. 

- :meth:`extract_region`: Extracts a local region of the surface, which is useful for performing operations on a small portion of the surface without affecting the full surface. This returns a `LocalSurface`` object.
- :meth:`add_data`: Adds a new data variable that is associated with either faces (default) or nodes of the surface. 
- :meth:`update_elevation`: Updates the elevation data of the surface using the `new_elevation` argument. The data is added to the existing elevation data by default, unless the user specifies `overwrite=True`.
- :meth:`apply_diffusion`: Applies diffusive degradation to the surface, which is a process that smooths out the surface by averaging the elevations of neighboring faces. This is useful for simulating the effects of erosion over time.
- :meth:`slope_collapse`: Applies diffusion only to surfaces that are steeper than a given threshold  given by the argument `critical_slope_angle`, which is set to 35 degrees by default (a typical value for the angle of repose). 
- :meth:`apply_noise`: Applies tubulence-style simplex noise to the surface. This can be useful for generating realistic surface features, such as hills and valleys. The noise is applied to the elevation data of the surface.
- :meth:`calculate_face_and_node_distances`: Returns one array for distances between a point and all faces, and another array with distances between that point and all nodes.
- :meth:`calculate_face_and_node_bearings`: Returns two arrays: Bearings (direction in degrees) between a point and all faces, Bearing between a point and all nodes.
- :meth:`find_nearest_face`: Returns the index of the face that is closest to a given point. 
- :meth:`find_nearest_node`: Returns the index of the node that is closest to a given point. 
- :meth:`interpolate_node_elevation_from_faces`: Some operations only affect faces, and so this method can be used to interpolate the elevation of nodes from the face elevations.
- :meth:`elevation_to_cartesian`: Convert elevation values to Cartesian coordinates. This is a basic utility function that takes the cartesian positions of the surface of a sphere and elevation values and returns the Cartesian coordinates of the surface with the elevations applied. This is used because the mesh is never altered in Cratermaker, but rather the elevations are applied to the mesh when it is visualized.
- :meth:`get_random_location_on_face`: Given a face index, this method will return a random location on that face. This is used to generate craters on a surface, as some surfaces have highly variable faces and therefore production functions must be evaluated on a face-by-face basis.

.. note::
    The LocalSurface object comes with distances and bearings from the center of the local region pre-computed. See below:

.. ipython:: python
    :okwarning:

    region=surface.extract_region(location=(205,45),region_radius=10e3)
    print(f"Region face distances:\n{region.face_distance}")
    print(f"Region face bearings:\n{region.face_bearing}")


Examples
--------

Extracting a local subset of the grid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose wish to extract a 10 km radius local region of the surface of the Moon centered at 45° N latitude, 205° E longitude:

.. ipython:: python
    :okwarning:

    surface=Surface.maker()
    region=surface.extract_region(location=(205,45), region_radius=10e3)
    print(f'Local region: {region}')

The ``region`` object now contains a view of all faces (along with their corresponding nodes and edges) of a local subset of the grid. Because it is a view of the surface not a copy, it allows for fast computation on small portions of the full grid..


Distances and bearings 
~~~~~~~~~~~~~~~~~~~~~~

The :class:`LocalSurface` object has two attributes called `face_distance` and `face_bearing` (and also `node_distance` and `node_bearing`), which are pre-computed distances and bearings from the center of the local region to all faces (nodes). This is useful for quickly accessing these values without having to compute them yourself.

.. ipython:: python
    :okwarning:

    print(f'Face distances from center of local region: {region.face_distance}')
    print(f'Face bearings from center of local region: {region.face_bearing}')

We can see that the local region contains points outside of the 10 km region. This is because  it contains a "buffer" of all faces that surround the outermost border of the local region, such that any operations that require neighboring faces across included edges or nodes can have access to them. We can see that this region contains 20 faces and 25 nodes, but only 7 of the faces actually lie within the 10 km local region. 
The other 13 are the buffer. From here, we can perform many of the calculations as seen in the list above. 

We can also calculate these for any arbitrary point within the local region (or full surface) using the methods :meth:`calculate_face_and_node_distances` and :meth:`calculate_face_and_node_bearings`. For instance, suppose we want to find the distances and bearings between a point at (205,45) and all faces and nodes in the local region:

.. ipython:: python
    :okwarning:

    face_distance, node_distance=region.calculate_face_and_node_distances(location=(205,45))
    print(f'Distances betwen location (205,45) and faces :{face_distance}')
    print(f'Distances betwen location (205,45) and nodes :{node_distance}')

With this method, two arrays are returned where the first array gives us an array of distances between the input location and the all faces, and the second array returns the distances from the input locaion and all of the nodes. It is best to use this method on smaller regions due to the size of the arrays if used on entire surface. We can do a similar calculation, but rather finding the distances between a location and the faces and nodes, we find the bearings: 

.. ipython:: python
    :okwarning:
    
    face_bearing, node_bearing =region.calculate_face_and_node_bearings(location=(205,45))
    print(f'Bearings betwen location (205,45) and faces :{face_bearing}')
    print(f'Bearings betwen location (205,45) and nodes :{node_bearing}')

As you can see from above, we recieve two arrays, which are the same sizes as the previous examples. However, they now tell us the "initial bearing" (the direction relative to due North) between a point and all faces and all nodes. 


Finding a face index
~~~~~~~~~~~~~~~~~~~~

Suppose we would like to find a face correspnding to a particulal location. We can do the following:

.. ipython:: python
    :okwarning:

    face_index=surface.find_nearest_face(location=(205,45))
    print(f'Nearest face to (205,40): {face_index}')

As seen above, we recieve an integer that gives us the index to the nearest face. One caveat is that this method will return the index of the face in which its center is the closest to the input location. Due to the shapes of the faces, this may or may not correspond to the face that contains the input location. However, it should correspond to at least one of the faces that borders the one containing the input location. You can view which faces these are using one of the built-in connectivity arrays. In this case, :attr:`face_face_connectivity`` contains the array of faces that are connected to a particular face:


.. ipython:: python
    :okwarning:

    print(f'Neighboring faces of {face_index} include: {surface.face_face_connectivity[face_index]}')

There are corresponding methods for finding the nearest node, as well as connectivity arrays for nodes, such as :attr:`node_face_connectivity` (faces associated with each node) and :attr:`face_node_connectivity` (nodes associated with each face).

.. note::
    Due to the variable number of nodes and edges associated with each face, there will sometimes be unused elements of the connectivity arrays. These are set to a large negative number, and so filtering out only indices greater than or equal to 0 will give you the valid indices. 


Converting elevation to Cartesian coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now lets say we are given an array of cartesian coordinates and an array of elevation, and we wish to convert the elevation values to Cartesian coordinates. We can do this by calling on :meth:`elevation_to_cartesian`: 

.. ipython:: python
    :okwarning:
    
    position = np.array([[0.0, 0.0, 1.0],
                        [0.707, 0.0, 0.707],
                        [0.0, 1.0, 0.0]])
    elevation = np.array([0.01, -0.02, 0.00])
    cartesian=surface.elevation_to_cartesian(position=position, elevation=elevation)
    print(cartesian)

From the result, we are returned an array the same size as the original poistion array, but with the elevation taken in account. Wherever there is a value other than 0 in the position array, the elevation is added to that value. 

More Surface examples
---------------------

See more examples at  :ref:`gal-surface`

