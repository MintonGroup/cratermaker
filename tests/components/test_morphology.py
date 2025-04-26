import unittest
import tempfile
from cratermaker.components.morphology import Morphology
from cratermaker.core.target import Target
from cratermaker.core.surface import Surface
from cratermaker import Crater

morphology_models = Morphology.available()
class TestMorphology(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.data_dir = self.temp_dir.name
        self.target = Target(name="Moon")
        self.pix = self.target.radius / 10.0
        self.gridlevel = 4
        self.surf = Surface.make(data_dir=self.data_dir, target=self.target, reset_surface=True, gridlevel=self.gridlevel)

        self.dummy_crater = Crater.make(
            location=(0.0, 0.0),
            final_diameter=100000.0,
        )

    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup()

    def test_model_registration(self):
        models = Morphology.available()
        self.assertIn("simplemoon", models)

    def test_model_instantiation(self):
        for model_name in morphology_models:
            morphology = Morphology.make(model_name)
            self.assertIsInstance(morphology, object)
            self.assertEqual(morphology.name, model_name)
            morphology.crater = self.dummy_crater
            self.assertIs(morphology.crater, self.dummy_crater)

    def test_invalid_crater_type_raises(self):
        with self.assertRaises(TypeError):
            for model_name in morphology_models:
                morphology = Morphology.make(model_name)
                morphology.crater = "not_a_crater"

    def test_form_crater_executes(self):
        for model_name in morphology_models:
            morphology = Morphology.make(model_name)
            morphology.form_crater(self.surf, crater=self.dummy_crater)

    def testmake_morphology(self):
        # Test the make_morphology function
        for model_name in morphology_models:
            morphology = Morphology.make(moprhology=model_name)
            self.assertIsInstance(morphology, Morphology)

if __name__ == '__main__':
    unittest.main()