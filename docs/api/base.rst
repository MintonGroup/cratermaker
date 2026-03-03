.. _api-base:

####
Base
####

The base classes define common methods and attributes that are shared across all components, such as the output directory where data files are stored. It also defines the interface for saving and loading data, which is used by all components to manage their data files.

Common Arguments
================

The `CommonArgs` class defines a set of common arguments that can be passed to any components. This is used to pass arguments to different components when they are created as part of a :py:class:`~cratermaker.core.simulation.Simulation` object.

.. autoclass:: cratermaker.core.base.CommonArgs
    :members:
    :undoc-members:
    :no-index-entry:

CratermakerBase
===============

This is a base class that defines common methods and attributes for all components in the cratermaker package. It includes methods for saving and loading data, as well as a reference to the output directory where data files are stored.

.. autoclass:: cratermaker.core.base.CratermakerBase
    :members:
    :undoc-members:
    :no-index-entry:

ComponentBase
==============

This is a base class for all components in the cratermaker package. It inherits from `CratermakerBase` and defines additional methods and attributes that are specific to high level component classes. 

.. autoclass:: cratermaker.core.base.ComponentBase
    :members:
    :undoc-members:
    :no-index-entry:


.. currentmodule:: cratermaker.core.base

.. autosummary::
   :toctree: generated/

    import_components

