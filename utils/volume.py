"""Render voxel cutaways of spatial simulations."""
import argparse
from pathlib import Path
import numpy as np
from utils.plot import plt
from matplotlib.ticker import MaxNLocator
from utils.results import load_run


def plot_volumes(directories, output):
    runs = [load_run(directory) for directory in directories]
    figure = plt.figure(figsize=(5 * len(runs), 6), layout="constrained")
    depth = max((np.nonzero(run["material"][run["meta"]["surface_z"]:] == 0)[0].max(initial=0) + 1)
                * run["meta"]["cell_nm"] for run in runs)
    for number, run in enumerate(runs, 1):
        m = run["meta"]
        if m.get("dimension", 2) != 3:
            raise ValueError("Voxel cutaways require spatial simulation output")
        start = m["surface_z"] - m["mask_cells"]
        stop = min(m["nz"], m["surface_z"] + int(depth / m["cell_nm"]) + 8)
        material = run["material"][start:stop].copy()
        material[:, :m["ny"] // 2, :] = 0
        occupied = (material != 0).transpose(2, 1, 0)
        colors = np.where(material.transpose(2, 1, 0) == 2, "#283a4b", "#dce8ee")
        x = (np.arange(m["nx"] + 1) - m["nx"] / 2) * m["cell_nm"]
        y = (np.arange(m["ny"] + 1) - m["ny"] / 2) * m["cell_nm"]
        z = (np.arange(start, stop + 1) - m["surface_z"]) * m["cell_nm"]
        xyz = np.meshgrid(x, y, z, indexing="ij")
        ax = figure.add_subplot(1, len(runs), number, projection="3d")
        ax.voxels(*xyz, occupied, facecolors=colors, edgecolor="none", shade=True)
        arguments = dict(line.split(" ", 1) for line in
                         (run["directory"] / "arguments.txt").read_text().splitlines())
        titles = {"cosine": r"$F_A=\cos\theta$", "enhanced": r"$F_B=(1+4\sin^2\theta)\cos\theta$",
                  "not_used": "Single-Cl chlorination"}
        ax.set(title=titles[arguments["--angular"]], xlabel="x (nm)", ylabel="y (nm)", zlabel="Depth (nm)")
        ax.set_box_aspect((m["nx"], m["ny"], stop - start))
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.zaxis.set_major_locator(MaxNLocator(5))
        ax.invert_zaxis()
        ax.view_init(elev=22, azim=-55)
    figure.suptitle("Spatial voxel geometry | front half removed for viewing")
    figure.savefig(output, dpi=160)
    plt.close(figure)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    plot_volumes(args.directories, args.output)
