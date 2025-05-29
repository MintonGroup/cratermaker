.. currentmodule:: cratermaker

.. image:: ../_static/surface_grid.png
    :alt: Surface
    :align: center
    :width: 300px

.. _ug-surface:

Surface
=======

Cratermaker's :ref:`Surface <api-surface>` component is used to represent target body's topography and other properties of its surface. Its prupose is to handle surface-related data by providing methods for setting elevation data, and calculating surface-related questions. This tool contains three classes of surface implementations that can be chosen by the user: **icosphere**, **arbitrary resolution**, **high resolution local**. 

 

Setting up a Surface object
---------------------------
All cratermaker components have a default configuration that is set when no arguments are passed through its 'maker' method. We can observe its defauts for Surfacing:

.. ipython:: python
   :okwarning:

   from cratermaker import Surface
   surface=Surface.maker()
   print(surface)

The default surface is automatically set to 'icosphere' with a grid level of 8 and the :ref:`Target <ug-target>` being the Moon. This surface consists of a uniform grid configuration with triangular faces that will be subdivided by the input value for the gridlevel argument.  Though it is limited to a few resolutions, it is the most accuracte and efficient class. As seen from the above, the effective pixel size is automatically set based on the surface and target that is used. The gridlevel arguement is useful because it controls the number of faces of the icosphere, and therefore defining its resolution via 10*4**gridlevel. We can see this change by doing the following:

.. ipython:: python
   :okwarning:

    surface_7=Surface.maker('icosphere',gridlevel=7)
    surface_9=Surface.maker('icosphere', gridlevel=9)
    print(f"Number of faces for a grid level of 7: {surface_7.n_face}",f"Number of faces for a grid level of 9: {surface_9.n_face}")
    print(f"The areas of each face for grid level of 7:{surface_7.face_areas}",f"The areas of each face for grid level of 9:{surface_9.face_areas}",)


 The user has the capability of changing the surface class as well as the target. This can be done by doing the following: 

.. ipython:: python
   :okwarning:

   surface=Surface.maker(surface='arbitrary_resolution', target='Mercury')
   print(surface)


As seen above, the arbitrary resolution grid also creates a uniform grid configurations but allows the user to define their pixel size. This class allows for any resolution that the user wants. 


.. ipython:: python
   :okwarning:

   surface=Surface.maker("hireslocal", pix=50, local_radius=1e3, local_location=(0,9), superdomain_scale_factor=1000)
   print(surface)

The last class for surfaces is the Hi-res local grid. There are additional arguments that the user can define based of their requirments. These are the pixel size, the local radius of the desired region, the local location in degrees, and the superdomain scale factor. The last argument allows the smallest craters that are modeled outside of the local region to be those whose ejecta could just reach the boundary. This class is useful for modeling local patches on the surface at a high resolution.


Using a Surface object
----------------------
Once you have surface object, you are now able to perform numerous surface-related computations. 


- :meth:`calculate_distance`: Takes a longitude and latitude pair and computes the great circle distance. 
- :meth:`calculate_bearing`: Takes a longitude and latitude pair and computes the intitial bearing from one point to another on the surface of the sphere.=
- :meth:`find_nearest_index`: Takes a longitude and latitude pair and calculates the Haversine Distance for each face of the grid. You will get a tuple that tells you the index of the face with the minimum distance.
- :meth:`calculate_face_and_node_distances`: Calculates the distances between nodes and faces at a given location
- :meth:`calculate_face_and_node_bearings`: Calculates the initial bearing between nodes and faces at a given location
- :meth:`interpolate_node_elevation_from_faces`: Computes the elevation for each node by taking the area-weighted average of the surrounding faces.


.. ipython:: python
    :okwarning:

    surface = Surface.maker()

    haversine=surface.calculate_distance(-1.457,0.732,2.738,-1.102)
    print(f"Haversine Distance: {haversine}")

    surface_2 = Surface.maker('arbitrary_resolution', target='Mars')
    surface_nearest_index=surface_2.find_nearest_index((21.37,124.82))
    print(f"Nearest index: {surface_nearest_index}")

    surface_3=Surface.maker("hireslocal", pix=50, local_radius=1e3, local_location=(0,9), superdomain_scale_factor=1000)