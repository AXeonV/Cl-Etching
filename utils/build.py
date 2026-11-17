"""Build the shared engine, dimensional backends, and their tests."""
from pathlib import Path
import argparse
import subprocess

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "build"
COMMON = list((ROOT / "src/common").glob("*.cpp"))
PLANAR = list((ROOT / "src/planar").glob("*.cpp"))
SPATIAL = list((ROOT / "src/spatial").glob("*.cpp"))


def build(target="simulation"):
    BUILD.mkdir(exist_ok=True)
    sources, name = {
        "simulation": ([*COMMON, *PLANAR, *SPATIAL, *ROOT.glob("src/app/*.cpp")], "etch_sim.exe"),
        "numerics": ([*COMMON, ROOT / "tests/common/test_numerics.cpp"], "etch_numerics_tests.exe"),
        "planar": ([*COMMON, *PLANAR, ROOT / "tests/planar/test_core.cpp"], "etch_tests.exe"),
        "spatial": ([*COMMON, *SPATIAL, ROOT / "tests/spatial/test_core.cpp"], "etch_spatial_tests.exe"),
        "tracking": ([*COMMON, *PLANAR, *SPATIAL, ROOT / "tests/test_tracking.cpp"], "etch_tracking_tests.exe"),
        "geometry": ([*COMMON, *PLANAR, ROOT / "tests/planar/test_geometry.cpp"], "etch_geometry_tests.exe"),
    }[target]
    executable = BUILD / name
    inputs = [*sources, *ROOT.glob("src/**/*.hpp")]
    if target != "simulation":
        inputs.extend(ROOT.glob("tests/**/*.hpp"))
    if not executable.exists() or any(p.stat().st_mtime > executable.stat().st_mtime for p in inputs):
        subprocess.run(["g++", "-std=c++17", "-O3", "-Wall", "-Wextra", "-I", str(ROOT / "src"),
                        *map(str, sources), "-o", str(executable)], check=True, cwd=ROOT)
    return executable


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", action="store_true")
    args = parser.parse_args()
    print(build())
    for target in ["numerics", "planar", "spatial", *(["geometry", "tracking"] if args.geometry else [])]:
        subprocess.run([str(build(target))], check=True)
