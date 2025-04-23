import unittest
import tempfile
from cratermaker.components.morphology import MorphologyModel, available_morphology_models, get_morphology_model, _init_morphology
from cratermaker.core.target import Target
from cratermaker.core.surface import Surface
from cratermaker import make_crater

morphology_models = available_morphology_models()
class TestMorphology(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.data_dir = self.temp_dir.name
        self.target = Target(name="Moon")
        self.pix = self.target.radius / 10.0
        self.gridlevel = 4
        self.surf = Surface.initialize(data_dir=self.data_dir, target=self.target, reset_surface=True, gridlevel=self.gridlevel)

        self.dummy_crater = make_crater(
            location=(0.0, 0.0),
            final_diameter=100000.0,
        )

    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup()

    def test_model_registration(self):
        models = available_morphology_models()
        self.assertIn("simplemoon", models)

    def test_model_instantiation(self):
        for model_name in morphology_models:
            morphology = get_morphology_model(model_name)()
            self.assertIsInstance(morphology, object)
            self.assertEqual(morphology.model, model_name)
            morphology.crater = self.dummy_crater
            self.assertIs(morphology.crater, self.dummy_crater)

    def test_invalid_crater_type_raises(self):
        with self.assertRaises(TypeError):
            for model_name in morphology_models:
                morphology = get_morphology_model(model_name)()
                morphology.crater = "not_a_crater"

    def test_form_crater_executes(self):
        for model_name in morphology_models:
            morphology = get_morphology_model(model_name)()
            morphology.form_crater(self.surf, crater=self.dummy_crater)

    def test_init_morphology(self):
        # Test the _init_morphology function
        for model_name in morphology_models:
            morphology = _init_morphology(moprhology=model_name)
            self.assertIsInstance(morphology, MorphologyModel)

if __name__ == '__main__':
    unittest.main()