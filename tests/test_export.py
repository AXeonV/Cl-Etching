"""Check material layout and physical coordinates in ParaView exports."""
import json
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
import numpy as np
from utils.export import export_structure

ROOT=Path(__file__).resolve().parents[1]


class StructureExportTests(unittest.TestCase):
    def test_planar_and_spatial_cells(self):
        for dimension in [2,3]:
            with self.subTest(dimension=dimension), tempfile.TemporaryDirectory(dir=ROOT/"build") as temporary:
                directory=Path(temporary)
                meta={"dimension":dimension,"nx":4,"nz":3,"surface_z":1,"cell_nm":1.0}
                if dimension==3:meta["ny"]=2
                shape=(3,2,4) if dimension==3 else (3,4)
                material=(np.arange(np.prod(shape))%3).astype("u1").reshape(shape)
                material.tofile(directory/"material.bin")
                (directory/"summary.json").write_text(json.dumps(meta),encoding="utf-8")
                tree=ET.parse(export_structure(directory))
                image=tree.find("ImageData")
                self.assertEqual(image.attrib["WholeExtent"],"0 4 0 2 0 3" if dimension==3 else "0 4 0 0 0 3")
                self.assertEqual(image.attrib["Origin"],"-2 -1 -1" if dimension==3 else "-2 0 -1")
                self.assertEqual(image.attrib["Spacing"],"1 1 1")
                arrays=tree.findall(".//DataArray")
                self.assertEqual(len(arrays),1)
                self.assertEqual(arrays[0].attrib,{"type":"UInt8","Name":"material","format":"ascii"})
                self.assertEqual(tree.findall(".//PointData/DataArray"),[])
                np.testing.assert_array_equal(np.fromstring(arrays[0].text,sep=" ",dtype="u1"),material.ravel())
