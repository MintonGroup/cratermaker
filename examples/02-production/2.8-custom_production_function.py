"""
Define a custom production function
=====================================

.. rubric:: By David Minton

This example demonstrates how to create a custom production function by defining a class derived from the base Production class. The only requirements are that the following two abstract methods are given concrete definitions:

.. code-block:: python

    @abstractmethod
    def _chronology(
        self,
        time_start: FloatLike | ArrayLike,
        time_end: FloatLike | ArrayLike,
        **kwargs: Any,
    ) -> FloatLike | ArrayLike: ...


    @abstractmethod
    def _csfd(self, diameter: FloatLike | ArrayLike, **kwargs: Any) -> FloatLike | ArrayLike: ...


These are called by the public :py:meth:`~cratermaker.components.production.Production.chronology` and :py:meth:`~cratermaker.components.production.Production.csfd` methods after they have validated the inputs, so your implementations do not need to do any input checking. In particular, ``time_end`` will never be ``None`` when ``_chronology`` is called.

You then register a name for your production function and then can apply it to a Simulation.

For this example, we will generate an outer solar system production function based on Case A of Zahnle et al. (2003) [#]_.


References
----------

.. [#] Zahnle, K., Schenk, P., Levison, H., Dones, L., 2003. Cratering rates in the outer Solar System. Icarus 163, 263-289. `doi: 10.1016/S0019-1035(03)00048-4 <https://doi.org/10.1016/S0019-1035(03)00048-4>`_
"""

import numpy as np

from cratermaker import Production, Simulation


@Production.register("zahnle_case_A")
class ZahnleCaseA(Production):
    def _csfd(self, diameter, **kwargs):
        # TODO: implement the Zahnle et al. (2003) Case A CSFD
        raise NotImplementedError

    def _chronology(self, time_start, time_end, **kwargs):
        # TODO: implement the Zahnle et al. (2003) Case A chronology function
        raise NotImplementedError

    @property
    def generator_type(self) -> str:
        """
        Get the generator type of the production function.

        Returns
        -------
        str
            The generator type of the production function.
        """
        return "projectile"

    @generator_type.setter
    def generator_type(self, value: str) -> None:
        """
        The generator type for this production function is determined by the version.
        """
        raise NotImplementedError("The generator type cannot be changed for this production function.")
