.. _api-utils:

.. currentmodule:: cratermaker.utils

#################
Utility functions
#################

.. _api-utils-montecarlo:

Monte Carlo
===========

.. autosummary::
   :toctree: generated/

    montecarlo_utils.get_random_location
    montecarlo_utils.get_random_impact_angle
    montecarlo_utils.get_random_velocity
    montecarlo_utils.get_random_size
    montecarlo_utils.bounded_norm
    montecarlo_utils.get_random_location_on_face
    montecarlo_utils.get_random_impact_direction
    montecarlo_utils.sample_pikefit
    montecarlo_utils.sample_logfit_heteroskedastic

.. _api-utils-general:

General utilities
=================

.. autosummary::
    :toctree: generated/

    general_utils.Parameter
    general_utils.normalize_coords
    general_utils.validate_and_normalize_location
    general_utils.format_large_units
    general_utils.cleanup
    general_utils.parameter
    general_utils.convert_deg2m_res


.. _api-utils-visualization:

Visualization utilities
=======================

.. autosummary::
    :toctree: generated/

    visualization_utils.CratermakerTheme
    visualization_utils.toggle_pyvista_actor
    visualization_utils.update_pyvista_help_message


