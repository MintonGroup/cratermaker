.. currentmodule:: cratermaker

.. image:: ../_static/surface_grid.png
    :alt: Surface
    :align: center
    :width: 300px

.. _ug-surface:

Surface
=======

Cratermaker's :ref:`Surface <api-surface>` component is used to represent target body's topography and other properties of its surface. Its prupose is to handle surface-related data by providing methods for setting elevation data, and calculating surface-related questions. This tool contains three classes of surface implementations that can be chosen by the user: **icosphere**, **arbitrary resolution**, **high resolution local**. SurfaceView is another usefool tool that exists in the Surface object, in which allows numerous operations and a view of the regional grid insteads of the entire mesh. 

 

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

    surface_7=Surface.maker('icosphere',gridlevel=7)
    surface_9=Surface.maker('icosphere', gridlevel=9)
    print(f"Number of faces for a grid level of 7: {surface_7.n_face}",f"Number of faces for a grid level of 9: {surface_9.n_face}")
    print(f"The areas of each face for grid level of 7:{surface_7.face_areas}",f"The areas of each face for grid level of 9:{surface_9.face_areas}",)


As expected, we observe more faces with a lower area when we increase the grid level. The user has the capability of changing the surface class as well as the target. This can be done by doing the following: 

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


- :meth:`calculate_face_and_node_distances`: Computes the distances from a given location to all faces and nodes.
- :meth:`calculate_face_and_node_bearings`: Computes the initial bearing from a given location to all faces and nodes.
- :meth:`find_nearest_index`: Returns the indices of the face and node that are the closest to a given point. 


.. ipython:: python
    :okwarning:

    surface = Surface.maker()

    haversine=surface.calculate_distance(-1.457,0.732,2.738,-1.102)
    print(f"Haversine Distance: {haversine}")

    surface_2 = Surface.maker('arbitrary_resolution', target='Mars')
    surface_nearest_index=surface_2.find_nearest_index((21.37,124.82))
    print(f"Nearest index: {surface_nearest_index}")

    surface_3=Surface.maker("hireslocal", pix=50, local_radius=1e3, local_location=(0,9), superdomain_scale_factor=1000)
    bearing=surface_3.calculate_face_and_node_distances((-1.457,0.732))
    print(f"Face and NodeDistances:{bearing}")