import unittest
import cratermaker
from cratermaker import Projectile
from cratermaker import general_utils as gu
import numpy as np

class TestProjectile(unittest.TestCase):
    
    def setUp():
       self.sim = cratermaker.Simulation() 

    def test_set_defaults(self):
       pass
    
if __name__ == '__main__':
    unittest.main()