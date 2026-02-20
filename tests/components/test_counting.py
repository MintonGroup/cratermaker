import tempfile
import unittest

from numpy.random import default_rng

from cratermaker import Counting, Crater, Morphology, Surface, Target

target = Target(name="Moon")
pix = target.radius / 10000.0
local_location = (0, 0)
local_radius = pix * 50
superdomain_scale_factor = 10000
surface_args = {
    "icosphere": {
        "gridlevel": 6,
    },
    "hireslocal": {
        "pix": pix,
        "local_location": local_location,
        "local_radius": local_radius,
        "superdomain_scale_factor": superdomain_scale_factor,
    },
}

crater = {
    "icosphere": {"diameter": 1000.0e3, "location": local_location},
    "hireslocal": {"diameter": 1000.0, "location": local_location},
}


class TestCounting(unittest.TestCase):
    def test_emplace(self):
        for surface_type, args in surface_args.items():
            with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
                surface = Surface.maker(surface_type, simdir=simdir, target=target, ask_overwrite=False, reset=True, **args)
                counting = Counting.maker(surface=surface)
                morphology = Morphology.maker(counting=counting, ejecta_truncation=4.0)
                morphology.emplace(**crater[surface_type])
                self.assertEqual(counting.n_observed, 1, msg=f"Failed for surface type {surface_type}")
                self.assertEqual(counting.n_emplaced, 1, msg=f"Failed for surface type {surface_type}")
                del surface
                del counting
                del morphology

    def test_cookie_cutting(self):
        for surface_type, args in surface_args.items():
            with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
                surface = Surface.maker(surface_type, simdir=simdir, target=target, ask_overwrite=False, reset=True, **args)
                counting = Counting.maker(surface=surface)
                morphology = Morphology.maker(counting=counting, ejecta_truncation=4.0)
                diameter = crater[surface_type]["diameter"]
                crater1 = Crater.maker(diameter=diameter, location=local_location)
                crater2 = Crater.maker(diameter=2 * diameter, location=local_location)
                morphology.emplace(crater1)
                self.assertEqual(counting.n_observed, 1, msg=f"Failed for surface type {surface_type}")
                self.assertEqual(counting.n_emplaced, 1, msg=f"Failed for surface type {surface_type}")
                morphology.emplace(crater2)
                self.assertEqual(counting.n_observed, 1, msg=f"Failed for surface type {surface_type}")
                self.assertEqual(counting.n_emplaced, 2, msg=f"Failed for surface type {surface_type}")
                del surface
                del counting
                del morphology

    def test_visibility_function(self):
        pass
        # local_radius = 2000.0
        # pix = 10.0
        # local_location = (0, 0)
        # h0 = 100.0
        # sigma = 200.0
        # kdiff = 20000.0
        # superdomain_scale_factor = 10000.0
        # with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
        #     surface = Surface.maker(
        #         "hireslocal",
        #         simdir=simdir,
        #         local_location=local_location,
        #         pix=pix,
        #         local_radius=local_radius,
        #         superdomain_scale_factor=superdomain_scale_factor,
        #         ask_overwrite=False,
        #     )
        #     counting = Counting.maker("minton2019", surface=surface)
        #     morphology = Morphology.maker(counting=counting, ejecta_truncation=4)
        #     morphology.emplace(diameter=1000.0)


if __name__ == "__main__":
    unittest.main()
