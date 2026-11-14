"""Plot one or more saved runs without changing or averaging their geometry."""
import argparse
import os
from pathlib import Path
os.environ["MPLCONFIGDIR"] = str(Path(__file__).resolve().parents[1] / "build/matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
from utils.results import load_run


def plot_runs(directories, output):
    runs = [load_run(directory) for directory in directories]
    half_width = 0
    max_depth = 0
    for run in runs:
        m = run["meta"]
        z, x = np.nonzero(run["material"][m["surface_z"]:] == 0)
        half_width = max(half_width, np.max(np.abs(x + 0.5 - m["nx"] / 2)) * m["cell_nm"])
        max_depth = max(max_depth, (np.max(z) + 1) * m["cell_nm"])
    figure, axes = plt.subplots(1, len(directories), figsize=(4.8 * len(directories), 7), squeeze=False,
                                constrained_layout=True)
    for ax, directory, run in zip(axes[0], directories, runs):
        m = run["meta"]
        arguments = dict(line.split(" ", 1) for line in (Path(directory) / "arguments.txt").read_text().splitlines())
        angular = arguments["--angular"]
        color = "#007f9b" if angular == "cosine" else "#de7438" if angular == "enhanced" else "#7e5794"
        extent = [-m["nx"] * m["cell_nm"] / 2, m["nx"] * m["cell_nm"] / 2,
                  (m["nz"] - m["surface_z"]) * m["cell_nm"], -m["surface_z"] * m["cell_nm"]]
        ax.imshow(run["material"], cmap=ListedColormap(["white", "#dce8ee", "#283a4b"]),
                  vmin=0, vmax=2, extent=extent, interpolation="nearest")
        segments = []
        for row in run["curves"]:
            x, z, tx, tz, nx, nz, a, b, c, _ = row
            u = np.linspace(-np.sqrt(2), np.sqrt(2), 41)
            v = a*u*u + b*u + c
            local = np.column_stack((tx*u + nx*v, tz*u + nz*v))
            points = local + [x, z]
            points[(np.abs(local) > 1).any(axis=1)] = np.nan
            points -= [m["nx"] / 2, m["surface_z"]]
            segments.append(points * m["cell_nm"])
        ax.add_collection(LineCollection(segments, colors=color, linewidths=0.65))
        titles = {"cosine": "Yield model\n" + r"$F_A(\theta)=\cos\theta$",
                  "enhanced": "Yield model\n" + r"$F_B(\theta)=(1+4\sin^2\theta)\cos\theta$",
                  "not_used": "Single-Cl chlorination\n" + r"$\mathrm{Si}\rightarrow\cdots\rightarrow\mathrm{SiCl_4(g)}$"}
        ax.set(title=titles[angular], xlabel="Lateral position (nm)", ylabel="Depth from initial Si surface (nm)")
        ax.set_xlim(-half_width - 35, half_width + 35)
        ax.set_ylim(max_depth + 50, -m["mask_cells"] * m["cell_nm"] - 25)
        ax.axhline(0, color="#93a3ac", linestyle="--", linewidth=0.6)
        ax.text(0.03, 0.025, f"{m['particles']:,} incident particles", transform=ax.transAxes, fontsize=9)
        ax.set_aspect("equal")
    figure.suptitle("Fixed-input rerun: material grid and local quadratic fits", fontsize=15)
    figure.legend(handles=[Patch(facecolor="#dce8ee", label="Si cells"),
                           Patch(facecolor="#283a4b", label="Fixed mask"),
                           Line2D([0], [0], color="#007f9b", label="Local fitted collision curves")],
                  loc="outside lower center", ncol=3, frameon=False, fontsize=9)
    figure.savefig(output, dpi=160)
    plt.close(figure)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    plot_runs(args.directories, args.output)
