"""Exercise dimensional dispatch and output contracts through the executable."""
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from utils.results import load_run, inventory_errors

ROOT = Path(__file__).resolve().parents[1]


class CommandLineTests(unittest.TestCase):
    def run_case(self, dimension, model, directory, explicit=True, options=()):
        command = [str(ROOT / "build/etch_sim.exe"), "--nx", "24", "--nz", "64", "--surface-z", "16",
                   "--mask-thickness", "4", "--opening", "12", "--particles", "3000", "--fwhm", "10",
                   "--energy-ev", "100", "--model", model, "--max-bounces", "500", "--output", str(directory)]
        if explicit:
            command += ["--dimension", str(dimension)]
        if dimension == 3:
            command += ["--ny", "12", "--opening-y", "8"]
        command += list(options)
        subprocess.run(command, check=True, capture_output=True)
        return load_run(directory)

    def test_both_dimensions_and_models(self):
        with tempfile.TemporaryDirectory(dir=ROOT / "build") as temporary:
            for dimension in [2, 3]:
                for model in ["yield", "chlorination"]:
                    with self.subTest(dimension=dimension, model=model):
                        run = self.run_case(dimension, model, Path(temporary) / f"{dimension}-{model}")
                        self.assertEqual(run["meta"]["dimension"], dimension)
                        self.assertEqual(run["material"].shape, (64, 12, 24) if dimension == 3 else (64, 24))
                        self.assertEqual(run["sites"].shape, (*run["material"].shape, 4))
                        self.assertEqual(inventory_errors(run)["silicon"], 0)
                        if model == "chlorination":
                            self.assertEqual(inventory_errors(run)["chlorine"], 0)
                        self.assertTrue(np.isfinite(run["curves"]).all())
                        self.assertEqual(run["curves"].shape[1], 19 if dimension == 3 else 10)

    def test_spatial_distribution_options(self):
        with tempfile.TemporaryDirectory(dir=ROOT / "build") as temporary:
            for angle in ["normal", "cosine", "gaussian"]:
                with self.subTest(angle=angle):
                    run = self.run_case(3, "yield", Path(temporary) / angle,
                                        options=("--angle", angle, "--energy", "gaussian", "--energy-sigma", "5"))
                    self.assertEqual(inventory_errors(run)["silicon"], 0)
                    self.assertTrue(np.isfinite(run["curves"]).all())

    def test_legacy_arguments_default_to_planar(self):
        with tempfile.TemporaryDirectory(dir=ROOT / "build") as temporary:
            old = self.run_case(2, "yield", Path(temporary) / "default", explicit=False)
            explicit = self.run_case(2, "yield", Path(temporary) / "explicit")
            np.testing.assert_array_equal(old["material"], explicit["material"])
            np.testing.assert_array_equal(old["sites"], explicit["sites"])


if __name__ == "__main__":
    unittest.main()
