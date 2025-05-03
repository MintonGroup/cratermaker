import unittest
import tempfile
from cratermaker.components.morphology import Morphology
from cratermaker.components.target import Target
from cratermaker.components.surface import Surface
from cratermaker import Crater, Simulation

morphology_models = Morphology.available()
class TestMorphology(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.target = Target(name="Moon")
        self.pix = self.target.radius / 10.0
        self.gridlevel = 5

        self.dummy_crater = Crater.maker(
            location=(0.0, 0.0),
            final_diameter=100000.0,
        )

    def test_model_registration(self):
        models = Morphology.available()
        self.assertIn("simplemoon", models)

    def test_model_instantiation(self):
        for model_name in morphology_models:
            morphology = Morphology.maker(model_name)
            self.assertIsInstance(morphology, object)
            self.assertEqual(morphology.name, model_name)
            morphology.crater = self.dummy_crater
            self.assertIs(morphology.crater, self.dummy_crater)

    def test_invalid_crater_type_raises(self):
        with self.assertRaises(TypeError):
            for model_name in morphology_models:
                morphology = Morphology.maker(model_name)
                morphology.crater = "not_a_crater"

    def test_form_crater_executes(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel)
            for model_name in morphology_models:
                morphology = Morphology.maker(model_name)
                morphology.form_crater(surface, crater=self.dummy_crater)

    def testmake_morphology(self):
        # Test the make_morphology function
        for model_name in morphology_models:
            morphology = Morphology.maker(moprhology=model_name)
            self.assertIsInstance(morphology, Morphology)

    def test_crater_depth(self):
        # Tests that the surface elevations are expected

        final_diameter_list = [100e3, 200e3, 500e3, 1000e3]
        delta_vals = [0.4, 0.3, 0.3, 0.2]
            
        gridargs = {
            "icosphere": {
                "gridlevel" : 6
                },
            "arbitrary_resolution": {
                "pix": 10e3,
                },
            "hireslocal": {
                "pix": 1.0e3,
                "local_location": (0, 0), 
                "local_radius": 100e3, 
                }
            }

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            for name, args in gridargs.items():
                sim = Simulation(simdir=simdir, surface=name, **args)
                for final_diameter, delta in zip(final_diameter_list, delta_vals):
                    sim.surface.reset()
                    # verify that the surface is flat
                    self.assertAlmostEqual(sim.surface.node_elevation.min(), 0.0, delta=1e0)
                    self.assertAlmostEqual(sim.surface.face_elevation.min(), 0.0, delta=1e0)
                    self.assertAlmostEqual(sim.surface.node_elevation.max(), 0.0, delta=1e0)
                    self.assertAlmostEqual(sim.surface.face_elevation.max(), 0.0, delta=1e0)
                    
                    sim.emplace_crater(final_diameter=final_diameter, location=(0, 0))

                    # Verify that the crater depth and rim heights are close to the expected values
                    self.assertAlmostEqual(-sim.surface.node_elevation.min() / sim.morphology.floordepth, 1.0, delta=delta)
                    self.assertAlmostEqual(-sim.surface.face_elevation.min() / sim.morphology.floordepth, 1.0, delta=delta)
                    self.assertAlmostEqual(sim.surface.node_elevation.max() / sim.morphology.rimheight, 1.0, delta=2*delta)
                    self.assertAlmostEqual(sim.surface.face_elevation.max() / sim.morphology.rimheight, 1.0, delta=2*delta)

if __name__ == '__main__':
    unittest.main()