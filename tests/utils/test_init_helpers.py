import unittest
import tempfile
from pathlib import Path
import os

class TestInitHelpers(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.simdir = self.temp_dir.name
        self.cwd = Path.cwd()

    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup() 
        os.chdir(self.cwd)  # Change back to the original working directory
        return       


if __name__ == '__main__':
    unittest.main()
