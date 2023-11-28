import unittest
from unittest.mock import patch, MagicMock
import cratermaker 

class TestSimulation(unittest.TestCase):

    def test_simulation_defaults(self):
        sim = cratermaker.Simulation()
        self.assertEqual(sim.target.name, "Moon")
        self.assertEqual(sim.target.material.name, "Soft Rock")

if __name__ == '__main__':
    unittest.main()
