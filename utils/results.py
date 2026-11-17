"""Read simulation arrays and check their inventories."""
import json
from pathlib import Path
import numpy as np


def load_run(directory):
    directory = Path(directory)
    meta = json.loads((directory / "summary.json").read_text())
    shape = ((meta["nz"], meta["ny"], meta["nx"]) if meta.get("dimension", 2) == 3
             else (meta["nz"], meta["nx"]))
    run = {
        "directory": directory, "meta": meta,
        "material": np.fromfile(directory / "material.bin", dtype="u1").reshape(shape),
        "sites": np.fromfile(directory / "sites.bin", dtype="<u2").reshape(*shape, 4),
        "primary": np.fromfile(directory / "primary_hits.bin", dtype="<u4").reshape(shape),
        "secondary": np.fromfile(directory / "secondary_hits.bin", dtype="<u4").reshape(shape),
        "curves": np.loadtxt(directory / "surface.csv", delimiter=",", skiprows=1, ndmin=2),
    }

    field = directory / "surface_nodes.bin"
    if field.exists():
        resolution = meta["surface_subdivisions"]
        field_shape = tuple(size * resolution + 1 for size in shape)
        run["surface_nodes"] = np.fromfile(field, dtype="<f8").reshape(*field_shape, 1 << len(shape))
    return run


def inventory_errors(run):
    meta = run["meta"]
    silicon = int(run["sites"].sum()) + meta["removed_si"] + meta.get("pruned_si",0) + meta.get("detached_si",0) - meta["initial_si"]
    chlorine = (int((run["sites"] * np.arange(4)).sum()) + 4 * meta["sicl4_products"]
                + meta.get("pruned_cl",0) + meta.get("detached_cl",0) - meta["captures"]) if meta["model"] == "chlorination" else None
    return {"silicon": silicon, "chlorine": chlorine}
