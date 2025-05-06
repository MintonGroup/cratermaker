.. currentmodule:: cratermaker

.. image:: ../_static/synthetic_complex.png
    :alt: Morphology
    :align: center
    :width: 300px

.. _ug-morphology:

Morphology
==========

The Morphology component is used to alter the topography of a :ref:`Surface <ug-surface>` object using a :ref:`Crater <ug-crater>` object. There is currently only one morphology model available ("simplemoon"), which is derived from CTEM. Like other components, a :ref:`Morphology <api-morphology>` model is generated using its :meth:`maker` method:

.. ipython:: python

    from cratermaker import Morphology
    morphology = Morphology.maker("simplemoon")
    print(morphology)

The morphology model can't do much without a crater. You can either pass a crater object in as an argument to the :meth:`maker` method, or you can set the crater after the fact. 

.. ipython:: python

    from cratermaker import Morphology, Crater
    morphology = Morphology.maker(crater=Crater.maker(final_diameter=50e3))
    print(morphology)

Every time you set a new crater, the morphology model will recompute its parameters.

.. ipython:: python

    morphology.crater = Crater.maker(final_diameter=1e3)
    print(morphology)
    morphology.crater = Crater.maker(final_diameter=15e3)
    print(morphology)

The primary purpose of a morphology model is to alter the topography of a surface, so typically you would not be using the morphology model directly. However, the built in morphology creation methods, such as :meth:`morphology.crater_shape` and :meth:`ejecta_shape` can be used in other context. See the example gallery for some examples of the the morphology model in action.


More Morphology examples
------------------------

See more complex usage examples in the gallery: :ref:`gal-morphology`

.. toctree::
   :maxdepth: 2
   :hidden:
