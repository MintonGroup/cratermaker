.. currentmodule:: cratermaker

.. image:: ../_static/surface_grid.png
    :alt: Surface
    :align: center
    :width: 300px

.. _ug-surface:

Surface
=======

Cratermaker's :ref:`Surface <api-surface>` component is used to represent target body's topography and other properties of its surface. Its prupose is to handle surface-related data by providing methods for setting elevation data, and calculating surface-related questions. This tool contains three classes of surface implementations that can be chosen by the user: **icosphere**, **arbitrary resolution**, **high resolution local**. LocalSurface is another usefool tool that exists in the Surface object, in which allows numerous operations and a view of the regional grid insteads of the entire mesh. 

 

Setting up a Surface object
---------------------------
All cratermaker components have a default configuration that is set when no arguments are passed through its 'maker' method. We can observe its defauts for Surfacing:

.. ipython:: python
   :okwarning:

   from cratermaker import Surface
   surface=Surface.maker()
   print(surface)

The default surface is automatically set to 'icosphere' with a grid level of 8 and the :ref:`Target <ug-target>` being the Moon. This surface consists of a uniform grid configuration with triangular faces that will be subdivided by the input value for the gridlevel argument.  Though it is limited to a few resolutions, it is the most accuracte and efficient class. As seen from the above, the effective pixel size is automatically set based on the surface and target that is used. The gridlevel arguement is useful because it controls the number of faces of the icosphere, and therefore defining its resolution via :math:`10\times4^{gridlevel}`. We can see this change by doing the following:

.. ipython:: python
   :okwarning:

    number_faces_4=10*4**4
    print(f"Number of faces for a grid level of 4: {number_faces_4}")

    surface_area_sphere_4=4*3.14*1738000**2
    average_area_4=surface_area_sphere_4/number_faces_4
    print(f"Average face area for grid level 4: {average_area_4}")


    number_faces_5=10*4**5
    print(f"Number of faces for a grid level of 5: {number_faces_5}")

    surface_area_sphere_5=4*3.14*1738000**2
    average_area_5=surface_area_sphere_5/number_faces_5
    print(f"Average face area for grid level 5: {average_area_5}")

As expected, we observe more faces with a lower area when we increase the grid level. The number of faces and the average face area can be determined by using the .n_face command and by taking the mean of .face_areas. Users also have the capability of changing the surface class as well as the target. This can be done by doing the following: 

.. ipython:: python
   :okwarning:

   surface=Surface.maker(surface='arbitrary_resolution', target='Mercury')
   print(surface)


As seen above, the arbitrary resolution grid also creates a uniform grid configurations but allows the user to define their pixel size. Unlike the icosphere, the faces for this class can vary in shape, making it less ideal. However, this class is useful when the user wants more flexibility in the resolution regime. 


.. ipython:: python
   :okwarning:

   surface=Surface.maker("hireslocal", pix=50, local_radius=1e3, local_location=(0,9), superdomain_scale_factor=1000)
   print(surface)

The last class for surfaces is the Hi-res local grid. There are additional arguments that the user can define based of their requirments. These are the pixel size, the local radius of the desired region, the local location in degrees, and the superdomain scale factor. The last argument allows the smallest craters that are modeled outside of the local region to be those whose ejecta could just reach the boundary. This class is useful for modeling local patches on the surface at a high resolution. Thus, if the user wants to simulate an impactor and study its ejecta behvior and its effect on local topography, this class would be useful due to its high resolution capabilities. 


Using a Surface object
----------------------
Once you have surface object, you are now able to perform numerous surface-related computations. 


- :meth:`calculate_face_and_node_distances`: Returns one array for distances between a point and all faces, and another array with distances between that point and all nodes.
- :meth:`calculate_face_and_node_bearings`: Returns two arrays: Bearings (direction in degrees) between a point and all faces, Bearing between a point and all nodes.
- :meth:`find_nearest_index`: Returns the indices of the face and node that are the closest to a given point. 


.. ipython:: python
    :okwarning:

    surface = Surface.maker()

    surface_2 = Surface.maker('arbitrary_resolution', target='Mars')
    surface_nearest_index=surface_2.find_nearest_index((21.37,124.82))
    print(f"Nearest index: {surface_nearest_index}")

    region=surface.extract_region(location=[0,3],region_radius=30000)
    face_bearing, node_bearing=region.calculate_face_and_node_bearings((1,2))
    print(face_bearing,node_bearing)

A SurfaceView is beneficial for the last example because you are able to a view a portion of the crater that is of interest. If you did the calculation without extracting a region, the user will recieve two large arrays that does not inform the user of any particular crater or portion of the surface. 