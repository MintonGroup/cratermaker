import unittest

from numpy.random import default_rng

from cratermaker import Counting, Crater, Morphology, Surface, Target

target = Target(name="Moon")
pix = target.radius / 10.0
local_location = (0, 0)
local_radius = pix * 2
superdomain_scale_factor = 100
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


class TestCounting(unittest.TestCase):
    def test_emplace(self):
        for surface_type, surface_args in gridargs.items():
            surface = Surface.maker(target=target, surface_type=surface_type, ask_overwrite=False, reset=True, **surface_args)
            counting = Counting.maker(surface=surface)
            morphology = Morphology.maker(counting=counting)
            morphology.emplace(diameter=1000.0e3, location=(0.0, 0.0))
            self.assertEqual(counting.n_observed, 1)
            self.assertEqual(counting.n_emplaced, 1)


if __name__ == "__main__":
    unittest.main()
