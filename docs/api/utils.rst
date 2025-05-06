.. _api-utils:

.. currentmodule:: cratermaker.utils

Utility functions
=================

.. _api-utils-montecarlo:

Monte Carlo
-----------

.. autosummary::
   :toctree: generated/

    montecarlo_utils.get_random_location
    montecarlo_utils.get_random_impact_angle
    montecarlo_utils.get_random_velocity
    montecarlo_utils.get_random_size
    montecarlo_utils.bounded_norm

.. _api-utils-general:

General utilities
-----------------

.. autosummary::
    :toctree: generated/

    general_utils.Parameter
    general_utils.validate_and_normalize_location
    general_utils.normalize_coords
    general_utils.R_to_CSFD

.. _api-utils-export:

Export functions
----------------

.. autosummary::
    :toctree: generated/

    export.to_vtk
    export.make_circle_file

.. _api-ComponentAPI:


Component API utilities
-----------------------

.. autoclass:: cratermaker.utils.component_utils.ComponentBase
   :members:
   :undoc-members:
   :no-index:

