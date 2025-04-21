import unittest
import tempfile
import os
from cratermaker.components.morphology import available_morphology_models, get_morphology_model
from cratermaker.core.target import Target
from cratermaker.core.surface import Surface
from cratermaker import Crater

class TestMorphology(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.simdir = self.temp_dir.name
        self.target = Target(name="Moon")
        self.pix = self.target.radius / 10.0
        self.gridlevel = 4
        self.surf = Surface.initialize(target=self.target, reset_surface=True, simdir=self.simdir, gridlevel=self.gridlevel)
        os.chdir(self.temp_dir.name)

        self.dummy_crater = Crater(
            location=(0.0, 0.0),
            diameter=1000.0,
        )

    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup()

    def test_model_registration(self):
        models = available_morphology_models()
        self.assertIn("simplemoon", models)

    def test_model_instantiation(self):
        cls = get_morphology_model("simplemoon")
        model = cls(crater=self.dummy_crater)
        self.assertIs(model.crater, self.dummy_crater)
        self.assertEqual(model.morphology_type, "simple")

    def test_invalid_crater_type_raises(self):
        cls = get_morphology_model("simplemoon")
        with self.assertRaises(TypeError):
            cls(crater="not_a_crater")

    def test_form_crater_executes(self):
        cls = get_morphology_model("simplemoon")
        morphology = cls(crater=self.dummy_crater)
        morphology.form_crater(self.surf)

if __name__ == '__main__':
    unittest.main()