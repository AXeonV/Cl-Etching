"""Evaluate the saved Hermite collision surface for plotting."""
import itertools
import numpy as np


def evaluate(nodes, coordinates, spacing):
    dimension = len(coordinates)
    cells = [np.minimum(np.floor(np.asarray(coordinate) / spacing).astype(int),
                        nodes.shape[-axis - 2] - 2)
             for axis, coordinate in enumerate(coordinates)]
    basis = []
    for coordinate, cell in zip(coordinates, cells):
        u = np.asarray(coordinate) / spacing - cell
        basis.append((2*u**3 - 3*u**2 + 1, -2*u**3 + 3*u**2,
                      u**3 - 2*u**2 + u, u**3 - u**2))
    value = 0.0
    for corner, mask in itertools.product(range(1 << dimension), repeat=2):
        weight = 1.0
        indices = []
        for axis in range(dimension):
            endpoint, derivative = (corner >> axis) & 1, (mask >> axis) & 1
            weight = weight * basis[axis][endpoint + 2*derivative]
            if derivative:
                weight = weight * spacing
            indices.append(cells[axis] + endpoint)
        value = value + weight * nodes[(*reversed(indices), mask)]
    return value


def section(run, slice_y=None):
    meta = run["meta"]
    resolution = meta["surface_subdivisions"]
    x = np.linspace(0, meta["nx"], meta["nx"] * resolution * 4 + 1)
    z = np.linspace(0, meta["nz"], meta["nz"] * resolution * 4 + 1)
    xx, zz = np.meshgrid(x, z)
    if meta["dimension"] == 2:
        coordinates = (xx, zz)
    else:
        y = meta["ny"] // 2 if slice_y is None else slice_y
        coordinates = (xx, y + 0.5, zz)
    values = evaluate(run["surface_nodes"], coordinates, 1.0 / resolution)
    return ((x - meta["nx"] / 2) * meta["cell_nm"],
            (z - meta["surface_z"]) * meta["cell_nm"], values)
