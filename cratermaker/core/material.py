    

        

@dataclass
class Material:
    """
    Represents the material properties relevant to the crater simulation.

    This class defines various physical properties of the material involved in the cratering process.
    

    Attributes
    ----------
    name : str
        The name of the material. If the material is matched to one that is present in the catalogue, the rest of the properties will be retrieved for it unless specified. If the name is not known from the catalogue, then all other properties must be supplied and in order to build a custom material.
    Ybar : float
        The strength of the material, typically defined in Pa. 
    other_properties : dict
        Other relevant properties of the material.

    Methods
    -------
    set_properties(name, **kwargs):
        Add a custom property to the material.

    """

    # Define all valid properties for the Target object
    name: str = None
    K1: float = None
    mu: float = None
    Ybar: float = None
    density: float = None 

    config_ignore = ['catalogue']  # Instance variables to ignore when saving to file
    def __post_init__(self):
        # Define some default crater scaling relationship terms (see Richardson 2009, Table 1, and Kraus et al. 2011 for Ice) 
        material_properties = [
            "name",       "K1",     "mu",   "Ybar",     "density" 
        ]
        material_values = [
            ("Water",     2.30,     0.55,   0.0,        1000.0),
            ("Sand",      0.24,     0.41,   0.0,        1750.0),
            ("Dry Soil",  0.24,     0.41,   0.18e6,     1500.0),
            ("Wet Soil",  0.20,     0.55,   1.14e6,     2000.0),
            ("Soft Rock", 0.20,     0.55,   7.60e6,     2250.0),
            ("Hard Rock", 0.20,     0.55,   18.0e6,     2500.0),
            ("Ice",       15.625,   0.48,   0.0,        900.0), 
        ]        
        
        self.catalogue = util._create_catalogue(material_properties, material_values)
        
        # Set properties for the Material object based on the catalogue value)
        if self.name:
            self.set_properties(catalogue=self.catalogue, key=self.name)
        else:
            raise ValueError('No material defined!')    
        
        return    
    
    def set_properties(self, **kwargs):
        """
        Set properties of the current object based on the provided keyword arguments.

        This function is a utility to update the properties of the current object. The actual implementation of the 
        property setting is handled by the `util._set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.

        Returns
        -------
        None
            The function does not return a value.
        """         
        util._set_properties(self,**kwargs)
        return
