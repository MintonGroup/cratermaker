import unittest
from cratermaker import Target, Scaling

class TestScale(unittest.TestCase):
   
    def test_scale_from_catalogue(self):
        # Test creating a scaling from the catalogue
        scaling = Scaling.maker(material="Soft Rock")
        self.assertEqual(scaling.K1, 0.20)
        self.assertEqual(scaling.mu, 0.55)
        self.assertEqual(scaling.Ybar, 7.60e6)
        self.assertEqual(scaling.target.density, 2250.0)
        self.assertEqual(scaling.projectile.density, 2250.0)
        self.assertIsInstance(scaling.K1, float)
        self.assertIsInstance(scaling.mu, float)
        self.assertIsInstance(scaling.Ybar, float)
        self.assertIsInstance(scaling.target.density, float)
        self.assertIsInstance(scaling.projectile.density, float)
        return

    def test_scale_override_catalogue(self):
        # Test overriding a property from the catalogue
        scaling = Scaling.maker(material="Soft Rock", Ybar=0.0)
        self.assertEqual(scaling.Ybar, 0.0)
        return

    def test_custom_scale(self):
        # Test creating a custom scaling
        scaling = Scaling.maker(material="Flubber", K1=3.8, mu=0.1, Ybar=1e7, density=2000.0)
        self.assertEqual(scaling.K1, 3.8)
        self.assertEqual(scaling.mu, 0.1)
        self.assertEqual(scaling.Ybar, 1e7)
        self.assertEqual(scaling.target.density, 2000.0)
        return

    def test_incomplete_custom_scale(self):
        # Test incomplete custom scaling creation
        with self.assertRaises(ValueError):
            Scaling.maker(material="Flubber", density=2000.0)
        with self.assertRaises(ValueError):
            Scaling.maker(material="Blorp", mu=0.1, Ybar=1e7, density=2000.0)
        return
    
    def test_scale_override_catalogue(self):
        # Test overriding a property from the catalogue
        target = Target(name="Mars")
        scaling = Scaling.maker(target=target,material="Sand")
        self.assertEqual(scaling.material, "Sand")    

if __name__ == '__main__':
    unittest.main()
