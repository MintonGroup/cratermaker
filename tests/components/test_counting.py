import unittest

from numpy.random import default_rng

from cratermaker import Crater, Counting, Target

gridlevel = 6
target = Target(name="Moon")
pix = target.radius / 10.0
local_location = (0, 0)
local_radius = pix * 2
superdomain_scale_factor = 100
gridargs = {
    "icosphere": {
        "gridlevel": gridlevel,
    },
    "hireslocal": {
        "pix": pix,
        "local_location": local_location,
        "local_radius": local_radius,
        "superdomain_scale_factor": superdomain_scale_factor,
    },
}

class TestCounting(unittest.TestCase):
    def test_tally(self):
        for 


if __name__ == "__main__":
    unittest.main()
