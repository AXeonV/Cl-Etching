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
                        self.assertEqual(run["meta"]["direction_components"], 3)
                        self.assertEqual(run["material"].shape, (64, 12, 24) if dimension == 3 else (64, 24))
                        self.assertEqual(run["sites"].shape, (*run["material"].shape, 4))
                        self.assertEqual(inventory_errors(run)["silicon"], 0)
                        if model == "chlorination":
                            self.assertEqual(inventory_errors(run)["chlorine"], 0)
                        self.assertTrue(np.isfinite(run["curves"]).all())
                        self.assertTrue(np.isfinite(run["surface_nodes"]).all())
                        self.assertEqual(run["curves"].shape[1], 19 if dimension == 3 else 10)

    def test_distribution_options(self):
        with tempfile.TemporaryDirectory(dir=ROOT / "build") as temporary:
            for dimension in [2, 3]:
                for angle in ["normal", "cosine", "gaussian"]:
                    with self.subTest(dimension=dimension, angle=angle):
                        run = self.run_case(dimension, "yield", Path(temporary) / f"{dimension}-{angle}",
                                            options=("--angle", angle, "--energy", "gaussian", "--energy-sigma", "5", "--surface-subdivisions", "2"))
                        self.assertEqual(inventory_errors(run)["silicon"], 0)
                        self.assertTrue(np.isfinite(run["curves"]).all())

    def test_saved_collision_field_matches_flat_plane(self):
        from utils.surface import evaluate
        with tempfile.TemporaryDirectory(dir=ROOT / "build") as temporary:
            for dimension in [2, 3]:
                run = self.run_case(dimension, "yield", Path(temporary) / str(dimension),
                                    options=("--particles", "0", "--opening", "24", "--opening-y", "12",
                                             "--surface-subdivisions", "2"))
                z = np.linspace(16.1, 16.9, 9)
                coordinates = (12.37, z) if dimension == 2 else (12.37, 6.25, z)
                np.testing.assert_allclose(evaluate(run["surface_nodes"], coordinates, 0.5),
                                           16.5 - z, atol=1e-11)

    def test_cell_size_parameter_preserves_grid_state(self):
        with tempfile.TemporaryDirectory(dir=ROOT / "build") as temporary:
            for dimension in [2,3]:
                with self.subTest(dimension=dimension):
                    original=self.run_case(dimension,"yield",Path(temporary)/f"{dimension}-coarse")
                    scaled=self.run_case(dimension,"yield",Path(temporary)/f"{dimension}-fine",
                                         options=("--cell-nm","1"))
                    self.assertEqual(original["meta"]["cell_nm"],2.5)
                    self.assertEqual(scaled["meta"]["cell_nm"],1.0)
                    np.testing.assert_array_equal(original["material"],scaled["material"])
                    np.testing.assert_array_equal(original["sites"],scaled["sites"])

    def test_inventory_update_invariants(self):
        with tempfile.TemporaryDirectory(dir=ROOT / "build") as temporary:
            for dimension in [2,3]:
                for model in ["yield","chlorination"]:
                    with self.subTest(dimension=dimension,model=model):
                        run=self.run_case(dimension,model,Path(temporary)/f"{dimension}-{model}")
                        self.assertEqual(run["meta"]["inventory_update"],"conservative")
                        self.assertEqual(run["meta"]["inventory_remap"],"slots")
                        self.assertEqual(run["meta"]["update_geometry"],"oriented_quadratic")
                        self.assertEqual(inventory_errors(run)["silicon"],0)
                        if model=="chlorination":self.assertEqual(inventory_errors(run)["chlorine"],0)
                        self.assertTrue((run["sites"][run["material"]==1].sum(axis=-1)>0).all())
                        self.assertTrue((run["sites"].sum(axis=-1)<=30).all())

    def test_default_dimension_is_planar(self):
        with tempfile.TemporaryDirectory(dir=ROOT / "build") as temporary:
            default = self.run_case(2, "yield", Path(temporary) / "default", explicit=False)
            explicit = self.run_case(2, "yield", Path(temporary) / "explicit")
            np.testing.assert_array_equal(default["material"], explicit["material"])
            np.testing.assert_array_equal(default["sites"], explicit["sites"])


if __name__ == "__main__":
    unittest.main()
