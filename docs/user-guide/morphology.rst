.. currentmodule:: cratermaker

.. image:: ../_static/synthetic_complex.png
    :alt: Morphology
    :align: center
    :width: 300px

.. _ug-morphology:

Morphology
==========

The Morphology component is used to alter the topography of a :class:`~cratermaker.Surface` object using a specified :class:`~cratermaker.Crater`. It supports multiple morphology models, such as "simplemoon", and is accessible via the class method :meth:`~cratermaker.Morphology.maker`.

.. ipython:: python

    from cratermaker.components.morphology import Morphology
    from cratermaker import Crater

    # Make a simple crater
    crater = Crater.maker(final_diameter=100e3)

    # Instantiate a morphology model (e.g., "simplemoon")
    morphology = Morphology.maker("simplemoon", crater=crater)

    # Access basic morphology outputs
    morphology.floor_depth, morphology.rim_height

You can also use Morphology without specifying a crater up front:

.. ipython:: python

    import numpy as np
    morphology = Morphology.maker("simplemoon")
    morphology.crater = Crater.maker(final_diameter=50e3)
    profile = morphology.crater_profile(np.array([0.0, 0.5, 1.0]))  # normalized radial profile



More Morphology examples
------------------------

See more complex usage examples in the gallery: :ref:`gal-morphology`

.. toctree::
   :maxdepth: 2
   :hidden:
