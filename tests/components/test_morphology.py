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

        final_diameter_list = [50e3, 100e3, 200e3, 500e3, 1000e3]
            
        gridargs = {
            "icosphere": {
                "level" :self.gridlevel, 
                },
            "arbitrary_resolution": {
                "pix": self.pix, 
                },
            "hireslocal": {
                "pix": self.pix, 
                "local_location": (0, 0), 
                "local_radius": 100e3, 
                "superdomain_scale_factor": 10
                }
            }

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            for name, args in gridargs.items():
                sim = Simulation(simdir=simdir, surface=name, **args)
                for final_diameter in final_diameter_list:
                    sim.emplace_crater(final_diameter=final_diameter, location=(0, 0))
                    node_depth = -float(sim.surface.node_elevation.min())
                    face_depth = -float(sim.surface.face_elevation.min())
                    self.assertAlmostEqual(node_depth, sim.morphology.floordepth, delta=1e2)
                    self.assertAlmostEqual(face_depth, sim.morphology.floordepth, delta=1e2)


if __name__ == '__main__':
    unittest.main()