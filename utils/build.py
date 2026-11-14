"""Build the simulator and its C++ unit tests."""
from pathlib import Path
import argparse
import subprocess

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "build"
CORE = [ROOT / "src" / f"{name}.cpp" for name in
        ("grid", "surface", "transport", "reaction", "cell_update", "simulation")]


def build(target="simulation"):
    BUILD.mkdir(exist_ok=True)
    main_path, executable_name = {
        "simulation": ("src/main.cpp", "etch_sim.exe"),
        "core": ("tests/test_core.cpp", "etch_tests.exe"),
        "geometry": ("tests/test_geometry.cpp", "etch_geometry_tests.exe"),
    }[target]
    main = ROOT / main_path
    executable = BUILD / executable_name
    sources = [*CORE, main]
    inputs = [*sources, *ROOT.glob("src/*.hpp")]
    if not executable.exists() or any(p.stat().st_mtime > executable.stat().st_mtime for p in inputs):
        subprocess.run(["g++", "-std=c++17", "-O3", "-Wall", "-Wextra", "-I", str(ROOT / "src"),
                        *map(str, sources), "-o", str(executable)], check=True, cwd=ROOT)
    return executable


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", action="store_true")
    args = parser.parse_args()
    print(build())
    subprocess.run([str(build("core"))], check=True)
    if args.geometry:
        subprocess.run([str(build("geometry"))], check=True)

