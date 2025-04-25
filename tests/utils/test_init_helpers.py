import unittest
from cratermaker.utils.init_helpers import _rng_init
import numpy as np
from numpy.random import Generator

class TestInitHelpers(unittest.TestCase):
    def test_argument_validation(self):

        # Invalid rng argument (not a Generator)
        with self.assertRaises(ValueError):
            _rng_init(rng=42)

        # Invalid rng_state type (not a dict)
        with self.assertRaises(ValueError):
            _rng_init(rng_state=123)

        # Invalid rng_seed (object instead of int)
        with self.assertRaises(ValueError):
            _rng_init(rng_seed="baz")

        # Accept random kwargs
        _rng_init(rng=np.random.default_rng(), foo="bar")

    def test_rng_consistency(self):
        # Basic initialization
        rng1, state1 = _rng_init()
        self.assertIsInstance(rng1, Generator)
        self.assertIsInstance(state1, dict)

        _ = rng1.uniform(0, 1)

        # Using the same rng instance
        rng2, _ = _rng_init(rng=rng1)
        _ = rng2.uniform(0, 1)

        # Using a fixed seed
        rng3, _ = _rng_init(rng_seed=42)
        u1 = rng3.uniform(0, 1)
        rng4, _ = _rng_init(rng_seed=42)
        u2 = rng4.uniform(0, 1)
        self.assertEqual(u1, u2)

        # Seed ignored when rng is provided
        rng5, _ = _rng_init(rng=rng3, rng_seed=123)
        u3 = rng5.uniform(0, 1)
        self.assertNotEqual(u3, u1)

        # Restore from state
        rng6, state1 = _rng_init(rng_seed=42)
        u4 = rng6.uniform(0, 1)
        rng7, _ = _rng_init(rng_state=state1)
        u5 = rng7.uniform(0, 1)
        self.assertEqual(u4, u5)

        # Modified state should diverge
        rng8, state2 = _rng_init(rng_seed=42)
        _ = rng8.uniform(0, 1)
        rng9, _ = _rng_init(rng=rng8, rng_seed=42)
        u6 = rng9.uniform(0, 1)
        self.assertNotEqual(u4, u6)


if __name__ == '__main__':
    unittest.main()
