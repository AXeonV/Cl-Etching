"""Export material cells as a ParaView ImageData structure."""
import argparse
import json
from pathlib import Path
import numpy as np


def export_structure(directory, output=None):
    directory = Path(directory)
    meta = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
    nx, nz = meta["nx"], meta["nz"]
    ny = meta["ny"] if meta["dimension"] == 3 else 0
    spacing = meta["cell_nm"]
    origin = (-nx * spacing / 2, -ny * spacing / 2, -meta["surface_z"] * spacing)
    material = np.fromfile(directory / "material.bin", dtype="u1").reshape(nz, max(ny, 1), nx)
    output = Path(output) if output is not None else directory / "structure.vti"
    extent = f"0 {nx} 0 {ny} 0 {nz}"
    xyz = " ".join(f"{value:.17g}" for value in origin)
    step = " ".join([f"{spacing:.17g}"] * 3)
    with output.open("w", encoding="utf-8", newline="\n") as file:
        file.write(f'''<?xml version="1.0"?>
<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">
  <ImageData WholeExtent="{extent}" Origin="{xyz}" Spacing="{step}">
    <Piece Extent="{extent}">
      <PointData/>
      <CellData Scalars="material">
        <DataArray type="UInt8" Name="material" format="ascii">
''')
        np.savetxt(file, material.reshape(-1, nx), fmt="%d")
        file.write("        </DataArray>\n      </CellData>\n    </Piece>\n  </ImageData>\n</VTKFile>\n")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    print(export_structure(args.directory, args.output))
