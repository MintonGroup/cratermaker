import tempfile
import unittest
from pathlib import Path

import numpy as np
from numpy.random import Generator

from cratermaker.core.base import (
    CommonArgs,
    CratermakerBase,
    _rng_init,
    _simdir_init,
    _to_config,
)


class TestBase(unittest.TestCase):
    def test_rng_argument_validation(self):
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

        # Using a fixed rng_seed
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
        rng8, _ = _rng_init(rng_seed=42)
        _ = rng8.uniform(0, 1)
        rng9, _ = _rng_init(rng=rng8, rng_seed=42)
        u6 = rng9.uniform(0, 1)
        self.assertNotEqual(u4, u6)

    def test_init_simdir(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            # Test with None
            newsimdir = _simdir_init(simdir)
            self.assertTrue(newsimdir.is_dir())

            # Test with a string path
            target_path = Path(simdir).resolve()
            newsimdir = _simdir_init(simdir=simdir)
            self.assertTrue(newsimdir.is_dir())
            self.assertEqual(newsimdir, target_path)

            # Test with a Path object
            simdir_path = _simdir_init(simdir=Path(simdir))
            self.assertTrue(simdir_path.is_dir())
            self.assertEqual(simdir_path, target_path)

            # Test with a relative path
            relative_path = "relative_simdir"
            newsimdir = _simdir_init(simdir=relative_path)
            self.assertTrue(newsimdir.is_dir())
            self.assertEqual(str(newsimdir), relative_path)
            newsimdir.rmdir()

            # Test with an invalid path
            with self.assertRaises(TypeError):
                _simdir_init(123)

    def test_to_config(self, **kwargs):
        class Dummy:
            def __init__(self):
                self._user_defined = [
                    "int_attr",
                    "float_attr",
                    "bool_attr",
                    "str_attr",
                    "np_int_attr",
                    "np_float_attr",
                    "np_array_attr",
                    "path_attr",
                    "dict_attr",
                    "list_attr",
                    "tuple_attr",
                    "none_attr",
                    "complex_attr",
                ]
                self.int_attr = 1
                self.float_attr = 2.5
                self.bool_attr = True
                self.str_attr = "hello"
                self.np_int_attr = np.int32(3)
                self.np_float_attr = np.float64(4.5)
                self.np_array_attr = np.array([1, 2, 3])
                self.path_attr = Path("/tmp/test")
                self.dict_attr = {"a": np.int16(10), "b": [np.float32(1.1), "x"]}
                self.list_attr = [1, np.float64(2.2)]
                self.tuple_attr = (3, np.int64(4))
                self.none_attr = None
                self.complex_attr = complex(1, 2)

        dummy = Dummy()
        config = _to_config(dummy)

        expected_keys = {
            "int_attr",
            "float_attr",
            "bool_attr",
            "str_attr",
            "np_int_attr",
            "np_float_attr",
            "np_array_attr",
            "path_attr",
            "dict_attr",
            "list_attr",
            "tuple_attr",
            "complex_attr",
        }

        self.assertEqual(set(config.keys()), expected_keys)
        self.assertEqual(config["int_attr"], 1)
        self.assertEqual(config["float_attr"], 2.5)
        self.assertEqual(config["bool_attr"], True)
        self.assertEqual(config["str_attr"], "hello")
        self.assertEqual(config["np_int_attr"], 3)
        self.assertEqual(config["np_float_attr"], 4.5)
        self.assertEqual(config["np_array_attr"], [1, 2, 3])
        self.assertEqual(config["path_attr"], str(Path("/tmp/test")))
        self.assertEqual(config["dict_attr"]["a"], 10)
        self.assertAlmostEqual(config["dict_attr"]["b"][0], 1.1, places=6)
        self.assertEqual(config["dict_attr"]["b"][1], "x")
        self.assertEqual(config["list_attr"][0], 1)
        self.assertAlmostEqual(config["list_attr"][1], 2.2, places=6)
        self.assertEqual(config["tuple_attr"], (3, 4))
        self.assertEqual(config["complex_attr"], str(complex(1, 2)))
        self.assertNotIn("none_attr", config)

    def test_cratermakerbase_and_commonargs(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            # Set up test values
            simdir = Path(simdir).resolve()
            rng, rng_state = _rng_init(rng_seed=123)

            # Test CommonArgs dataclass
            args = CommonArgs(simdir=simdir, rng=rng, rng_seed=123, rng_state=rng_state)
            self.assertEqual(args.simdir, simdir)
            self.assertEqual(args.rng_seed, 123)
            self.assertEqual(args.rng, rng)
            self.assertEqual(args.rng_state, rng_state)

            # Test CratermakerBase initialization
            base = CratermakerBase(
                simdir=simdir, rng=rng, rng_seed=123, rng_state=rng_state
            )

            # Attributes
            self.assertEqual(base.simdir, simdir)
            self.assertEqual(base._rng_seed, 123)
            self.assertEqual(base._rng, rng)  # direct access
            self.assertEqual(base._rng_state, rng_state)
            self.assertIsInstance(base.rng, Generator)
            self.assertEqual(base.rng_state, rng_state)

            # Test RNG reproducibility
            val1 = base.rng.uniform()
            new_rng, _ = _rng_init(rng_state=rng_state)
            val2 = new_rng.uniform()
            self.assertAlmostEqual(val1, val2, places=7)


if __name__ == "__main__":
    unittest.main()
