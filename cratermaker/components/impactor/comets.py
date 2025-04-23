import numpy as np
from numpy.typing import NDArray
from typing import Any
from cratermaker.utils import montecarlo as mc
from cratermaker.components.impactor import ImpactorModel, register_impactor_model

@register_impactor_model("comets")
class CometImpactors(ImpactorModel):
    """
    An operations class for computing the impactor properties of an asteroid source population.

    Parameters
    ----------
    **kwargs : Any
        Additional keyword arguments to be passed to internal functions.
    """
    
    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)

