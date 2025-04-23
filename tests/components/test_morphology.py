import unittest
import tempfile
from cratermaker.components.morphology import available_morphology_models, get_morphology_model
from cratermaker.core.target import Target
from cratermaker.core.surface import Surface
from cratermaker import Crater

class TestMorphology(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.data_dir = self.temp_dir.name
        self.target = Target(name="Moon")
        self.pix = self.target.radius / 10.0
        self.gridlevel = 4
        self.surf = Surface.initialize(data_dir=self.data_dir, target=self.target, reset_surface=True, gridlevel=self.gridlevel)

        self.dummy_crater = Crater(
            location=(0.0, 0.0),
            final_diameter=100000.0,
        )
        self.morphology = get_morphology_model("simplemoon")()

    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup()

    def test_model_registration(self):
        models = available_morphology_models()
        self.assertIn("simplemoon", models)

    def test_model_instantiation(self):
        self.morphology.crater = self.dummy_crater
        self.assertIs(self.morphology.crater, self.dummy_crater)

    def test_invalid_crater_type_raises(self):
        with self.assertRaises(TypeError):
            self.morphology.crater = "not_a_crater"

    def test_form_crater_executes(self):
        self.morphology.form_crater(self.surf, crater=self.dummy_crater)

if __name__ == '__main__':
    unittest.main()