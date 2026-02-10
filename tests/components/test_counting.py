import tempfile
import unittest

from numpy.random import default_rng

from cratermaker import Counting, Crater, Morphology, Surface, Target

target = Target(name="Moon")
pix = target.radius / 10000.0
local_location = (0, 0)
local_radius = pix * 50
superdomain_scale_factor = 10000
gridargs = {
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
        for surface_type, surface_args in gridargs.items():
            with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
                surface = Surface.maker(surface_type, simdir=simdir, target=target, ask_overwrite=False, reset=True, **surface_args)
                counting = Counting.maker(surface=surface)
                morphology = Morphology.maker(counting=counting, ejecta_truncation=2.0)
                morphology.emplace(**crater[surface_type])
                self.assertEqual(counting.n_observed, 1, msg=f"Failed for surface type {surface_type}")
                self.assertEqual(counting.n_emplaced, 1, msg=f"Failed for surface type {surface_type}")
                del surface
                del counting
                del morphology


if __name__ == "__main__":
    unittest.main()
