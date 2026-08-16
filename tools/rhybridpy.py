#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright 2013-2026 Finnish Meteorological Institute
# Copyright 2017-2026 University of Helsinki
#
# This single-file adaptation is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed without any warranty; without even the implied
# warranty of merchantability or fitness for a particular purpose.  See
# <https://www.gnu.org/licenses/> for details.
#
"""RHybridPy is based on the Analysator package (https://doi.org/10.5281/zenodo.4462514).

RHybridPy: a Python package for RHybrid VLSV workflows.

Usage:

    import rhybridpy as rhb

    vr = rhb.vlsvfile.VlsvReader("bulk.vlsv")
    values = vr.read_variable("rho")
    xyz, cellids, data, header = rhb.calculations.vlsv_intpol_points(
        vr, [[0.0, 0.0, 0.0]], ["rho"], interpolation_order=1
    )

Only NumPy and the Python standard library are required.

VLSV bulk arrays are stored in file order.  In parallel files that order is
the concatenation of the MPI ranks' local arrays, so load balancing may change
it between files.  ``CellID`` is therefore the stable key: pass ``cellids=``
when values from different files must be aligned.  Geometry and interpolation
methods perform this CellID-based alignment internally.
"""

from __future__ import annotations

import ast
import os
import shutil
import struct
import types
import warnings
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Any, Iterable, Sequence

import numpy as np

__all__ = [
    "VlsvReader",
    "VlsvParticles",
    "VlsvWriter",
    "VariableInfo",
    "vlsv_intpol_points",
    "vlsvfile",
    "calculations",
]
__version__ = "1.1.4"


def _literal(value: str | None, default: Any = None) -> Any:
    if value is None:
        return default
    return ast.literal_eval(value)


def _dtype(datatype: str, datasize: int) -> np.dtype:
    table = {
        ("float", 4): np.dtype("<f4"),
        ("float", 8): np.dtype("<f8"),
        ("int", 4): np.dtype("<i4"),
        ("int", 8): np.dtype("<i8"),
        ("uint", 4): np.dtype("<u4"),
        ("uint", 8): np.dtype("<u8"),
    }
    try:
        return table[(datatype.lower(), int(datasize))]
    except KeyError as exc:
        raise TypeError(
            f"Unsupported VLSV datatype/datasize: {datatype!r}/{datasize!r}"
        ) from exc


def _apply_operator(
    data: Any, operator: str | int | None, vector_size: int | None = None
) -> Any:
    array = np.asarray(data)
    if operator is None or operator == "pass":
        return data
    if isinstance(operator, (int, np.integer)):
        return array[..., int(operator)]
    op = str(operator).lower()
    if op in {"x", "0"}:
        return array[..., 0]
    if op in {"y", "1"}:
        return array[..., 1]
    if op in {"z", "2"}:
        return array[..., 2]
    if op == "magnitude":
        return (
            np.abs(array)
            if vector_size == 1 or array.ndim == 0
            else np.linalg.norm(array, axis=-1)
        )
    if op == "absolute":
        return np.abs(array)
    raise ValueError(f"Unsupported operator {operator!r}")


def _read_xml_footer(filename: str) -> ET.Element:
    with open(filename, "rb") as stream:
        stream.seek(8)
        raw = stream.read(8)
        if len(raw) != 8:
            raise IOError(f"{filename!r} is too short to be a VLSV file")
        (offset,) = struct.unpack("=Q", raw)
        stream.seek(offset)
        xml_data = stream.read()
    try:
        return ET.fromstring(xml_data)
    except ET.ParseError as exc:
        raise IOError(f"Invalid VLSV XML footer in {filename!r}") from exc


@dataclass(init=False)
class VariableInfo:
    """Data and metadata associated with a VLSV variable."""

    data: Any
    name: str
    units: str
    latex: str
    latexunits: str

    def __init__(
        self,
        data_array: Any = None,
        name: str = "",
        units: str = "",
        latex: str = "",
        latexunits: str = "",
        *,
        data: Any = None,
    ) -> None:
        if data_array is not None and data is not None:
            raise TypeError("Pass either data_array or data, not both")
        self.data = data if data is not None else data_array
        self.name = name
        self.units = units
        self.latex = latex or name.replace("_", r"\_")
        self.latexunits = latexunits


class VlsvReader:
    """Read stored variables and parameters from a uniform-grid VLSV file."""

    def __init__(self, file_name: str, **_: Any) -> None:
        self.file_name = os.path.abspath(os.fspath(file_name))
        if not os.path.isfile(self.file_name):
            raise FileNotFoundError(self.file_name)
        self.__xml_root = _read_xml_footer(self.file_name)
        self._cell_indices: dict[str, dict[int, int]] = {}
        self._cellid_arrays: dict[str, np.ndarray] = {}
        self._cellid_base: int | None = None
        self._parameters: dict[str, Any] = {}
        self._read_spatial_geometry()

    def __enter__(self) -> "VlsvReader":
        return self

    def __exit__(self, *_: Any) -> None:
        return None

    @property
    def xml_root(self) -> ET.Element:
        """The parsed footer, exposed read-only for compatible writers."""
        return self.__xml_root

    def _elements(
        self, *, name: str = "", tag: str = "", mesh: str | None = ""
    ) -> Iterable[ET.Element]:
        wanted_name = name.lower()
        matches: list[ET.Element] = []
        for child in self.__xml_root:
            if tag and child.tag != tag:
                continue
            if wanted_name:
                if "name" not in child.attrib:
                    continue
                if child.attrib["name"].lower() != wanted_name:
                    continue
            matches.append(child)

        if not mesh:
            yield from matches
            return

        # The original package requires an exact mesh match.  RHybridPy also
        # supports older meshless VLSV entries, but only as a fallback when no
        # exact entry exists.  This avoids mixing two differently permuted
        # CellID arrays when both legacy and mesh-tagged entries are present.
        exact = [child for child in matches if child.attrib.get("mesh") == mesh]
        if exact:
            yield from exact
        else:
            yield from [child for child in matches if "mesh" not in child.attrib]

    def _find(
        self, *, name: str = "", tag: str = "", mesh: str | None = ""
    ) -> ET.Element:
        try:
            return next(self._elements(name=name, tag=tag, mesh=mesh))
        except StopIteration as exc:
            raise KeyError(
                f"VLSV entry not found: name={name!r}, tag={tag!r}, mesh={mesh!r}"
            ) from exc

    def _read_spatial_geometry(self) -> None:
        try:
            bbox = np.asarray(
                self.read(tag="MESH_BBOX", mesh="SpatialGrid"), dtype=np.int64
            )
        except (KeyError, ValueError):
            bbox = None

        if bbox is not None and bbox.size >= 6:
            self._mesh_size = bbox[:3]
            self._block_size = bbox[3:6]
            x = np.asarray(self.read(tag="MESH_NODE_CRDS_X", mesh="SpatialGrid"))
            y = np.asarray(self.read(tag="MESH_NODE_CRDS_Y", mesh="SpatialGrid"))
            z = np.asarray(self.read(tag="MESH_NODE_CRDS_Z", mesh="SpatialGrid"))
            self._extent = np.array([x[0], y[0], z[0], x[-1], y[-1], z[-1]])
            return

        try:
            self._mesh_size = np.array(
                [
                    int(self.read_parameter("xcells_ini")),
                    int(self.read_parameter("ycells_ini")),
                    int(self.read_parameter("zcells_ini")),
                ],
                dtype=np.int64,
            )
            self._block_size = np.ones(3, dtype=np.int64)
            self._extent = np.array(
                [
                    self.read_parameter("xmin"),
                    self.read_parameter("ymin"),
                    self.read_parameter("zmin"),
                    self.read_parameter("xmax"),
                    self.read_parameter("ymax"),
                    self.read_parameter("zmax"),
                ],
                dtype=float,
            )
        except (KeyError, TypeError, ValueError):
            self._mesh_size = np.ones(3, dtype=np.int64)
            self._block_size = np.ones(3, dtype=np.int64)
            self._extent = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])

    def _read_element(self, child: ET.Element) -> np.ndarray | np.generic:
        vector_size = int(_literal(child.attrib.get("vectorsize"), 1))
        array_size = int(_literal(child.attrib.get("arraysize"), 1))
        element_size = int(_literal(child.attrib.get("datasize")))
        datatype = child.attrib["datatype"]
        offset = int(_literal(child.text))
        with open(self.file_name, "rb") as stream:
            stream.seek(offset)
            data = np.fromfile(
                stream,
                dtype=_dtype(datatype, element_size),
                count=array_size * vector_size,
            )
        if data.size != array_size * vector_size:
            raise IOError(
                f"Short VLSV read for {child.tag}/{child.attrib.get('name', '')}"
            )
        if vector_size > 1:
            return data.reshape(array_size, vector_size)
        if array_size == 1:
            return data[0]
        return data

    def read(
        self,
        name: str = "",
        tag: str = "",
        mesh: str | None = "",
        operator: str | int = "pass",
        cellids: Any = -1,
    ) -> Any:
        if not name and not tag:
            raise ValueError("read() requires a name or tag")
        child = self._find(name=name, tag=tag, mesh=mesh)
        data = self._read_element(child)
        if tag == "VARIABLE" and not _is_all_cells(cellids):
            data = self._select_cellids(data, cellids, child)
        vector_size = int(_literal(child.attrib.get("vectorsize"), 1))
        return _apply_operator(data, operator, vector_size)

    @staticmethod
    def _mesh_key(mesh: str | None) -> str:
        return mesh if mesh else "<legacy>"

    def _cellid_child_for(self, variable_child: ET.Element | None = None) -> ET.Element:
        mesh = variable_child.attrib.get("mesh") if variable_child is not None else "SpatialGrid"
        if mesh:
            return self._find(name="CellID", tag="VARIABLE", mesh=mesh)

        legacy = [
            child
            for child in self._elements(name="CellID", tag="VARIABLE")
            if "mesh" not in child.attrib
        ]
        if legacy:
            return legacy[0]
        return self._find(name="CellID", tag="VARIABLE", mesh="SpatialGrid")

    def _cellids_for(self, variable_child: ET.Element | None = None) -> np.ndarray:
        child = self._cellid_child_for(variable_child)
        key = self._mesh_key(child.attrib.get("mesh"))
        if key not in self._cellid_arrays:
            self._cellid_arrays[key] = np.atleast_1d(self._read_element(child)).astype(
                np.int64, copy=False
            )
        return self._cellid_arrays[key]

    def _ensure_cell_index(
        self, variable_child: ET.Element | None = None
    ) -> dict[int, int]:
        child = self._cellid_child_for(variable_child)
        key = self._mesh_key(child.attrib.get("mesh"))
        if key not in self._cell_indices:
            self._cell_indices[key] = {
                int(cellid): index
                for index, cellid in enumerate(self._cellids_for(variable_child))
            }
        return self._cell_indices[key]

    def _select_cellids(
        self, data: Any, cellids: Any, variable_child: ET.Element
    ) -> Any:
        index = self._ensure_cell_index(variable_child)
        scalar = np.isscalar(cellids)
        requested = np.atleast_1d(cellids).astype(np.int64)
        positions = []
        for cellid in requested:
            try:
                positions.append(index[int(cellid)])
            except KeyError as exc:
                raise KeyError(f"CellID {int(cellid)} is not stored in the file") from exc
        selected = np.asarray(data)[positions]
        return selected[0] if scalar else selected

    def get_all_variables(self) -> list[str]:
        return [
            child.attrib["name"]
            for child in self.__xml_root
            if child.tag == "VARIABLE" and "name" in child.attrib
        ]

    def get_variables(self) -> list[str]:
        return sorted(set(self.get_all_variables()))

    def check_variable(self, name: str) -> bool:
        wanted = name.lower()
        return any(variable.lower() == wanted for variable in self.get_all_variables())

    def check_parameter(self, name: str) -> bool:
        wanted = name.lower()
        return any(
            child.tag == "PARAMETER"
            and child.attrib.get("name", "").lower() == wanted
            for child in self.__xml_root
        )

    def list(
        self,
        parameter: bool = True,
        variable: bool = True,
        mesh: bool = False,
        datareducer: bool = False,
        operator: bool = False,
        other: bool = False,
    ) -> None:
        """Print a compact inventory of the stored VLSV content.

        ``datareducer`` is accepted for source compatibility but produces no
        entries because RHybridPy deliberately has no reducer registry.
        """
        if parameter:
            print("tag = PARAMETER")
            for child in self.__xml_root:
                if child.tag == "PARAMETER" and "name" in child.attrib:
                    print("   ", child.attrib["name"])
        if variable:
            print("tag = VARIABLE")
            for name in self.get_all_variables():
                print("   ", name)
        if mesh:
            print("tag = MESH")
            for child in self.__xml_root:
                if child.tag == "MESH" and "name" in child.attrib:
                    print("   ", child.attrib["name"])
        if datareducer:
            print("Datareducers: not included in RHybridPy")
        if operator:
            print("Data operators: pass, x, y, z, 0, 1, 2, magnitude, absolute")
        if other:
            print("Other:")
            for child in self.__xml_root:
                if child.tag not in {"PARAMETER", "VARIABLE", "MESH"}:
                    print(
                        "    tag =",
                        child.tag,
                        "mesh =",
                        child.attrib.get("mesh", ""),
                    )

    def read_attribute(
        self,
        name: str = "",
        mesh: str | None = "",
        attribute: str = "",
        tag: str = "",
    ) -> Any:
        """Return an XML-footer attribute for a matching VLSV entry."""
        if not name and not tag:
            raise ValueError("read_attribute() requires a name or tag")
        if not attribute:
            raise ValueError("read_attribute() requires an attribute name")
        child = self._find(name=name, tag=tag, mesh=mesh)
        if attribute not in child.attrib:
            raise KeyError(
                f"Attribute {attribute!r} is absent from "
                f"{tag or child.tag}/{name or child.attrib.get('name', '')}"
            )
        raw = child.attrib[attribute]
        try:
            return ast.literal_eval(raw)
        except (ValueError, SyntaxError):
            return raw

    def read_variable_vectorsize(self, name: str) -> int:
        return int(
            self.read_attribute(
                name=name,
                tag="VARIABLE",
                mesh="SpatialGrid",
                attribute="vectorsize",
            )
        )

    def read_metadata(
        self, name: str = "", tag: str = "VARIABLE", mesh: str | None = ""
    ) -> tuple[str, str, str, str]:
        child = self._find(name=name, tag=tag, mesh=mesh)
        return (
            child.attrib.get("unit", ""),
            child.attrib.get("unitLaTeX", ""),
            child.attrib.get("variableLaTeX", ""),
            child.attrib.get("unitConversion", ""),
        )

    def read_variable(
        self, name: str, cellids: Any = -1, operator: str | int = "pass"
    ) -> Any:
        """Read a SpatialGrid variable.

        With the default ``cellids=-1``, values are returned in this file's
        storage order.  That order can differ between VLSV files after load
        balancing.  Supplying CellIDs returns values in the requested CellID
        order and is the safe form for cross-file comparisons.
        """
        if not self.check_variable(name):
            raise KeyError(f"Variable {name!r} is not stored in {self.file_name!r}")
        child = self._find(name=name, tag="VARIABLE", mesh="SpatialGrid")
        mesh = child.attrib.get("mesh", "")
        return self.read(
            name=name,
            tag="VARIABLE",
            mesh=mesh,
            operator=operator,
            cellids=cellids,
        )

    def read_variable_info(
        self, name: str, cellids: Any = -1, operator: str | int = "pass"
    ) -> VariableInfo:
        child = self._find(name=name, tag="VARIABLE", mesh="SpatialGrid")
        units = child.attrib.get("unit", "")
        latexunits = child.attrib.get("unitLaTeX", "")
        latex = child.attrib.get("variableLaTeX", "")
        out_name = name if operator == "pass" else f"{name}_{operator}"
        if not latex:
            latex = name.replace("_", r"\_")
        if operator != "pass":
            latex = f"|{latex}|" if operator == "magnitude" else f"{latex}_{operator}"
        return VariableInfo(
            self.read_variable(name, cellids=cellids, operator=operator),
            name=out_name,
            units=units,
            latex=latex,
            latexunits=latexunits,
        )

    def read_parameter(self, name: str) -> Any:
        key = name.lower()
        if key in self._parameters:
            return self._parameters[key]
        aliases = ("t", "time") if key in {"t", "time"} else (name,)
        for candidate in aliases:
            try:
                value = self.read(name=candidate, tag="PARAMETER")
                self._parameters[key] = value
                return value
            except KeyError:
                pass
        raise KeyError(f"Parameter {name!r} is not stored in {self.file_name!r}")

    def get_spatial_mesh_extent(self) -> np.ndarray:
        return self._extent.copy()

    def get_spatial_mesh_size(self) -> np.ndarray:
        return self._mesh_size.copy()

    def get_spatial_block_size(self) -> np.ndarray:
        return self._block_size.copy()

    def get_cellid_locations(self) -> dict[int, int]:
        """Map each stored SpatialGrid CellID to its current file-array index.

        The indices are deliberately not assumed to be stable across files;
        they reflect the MPI-rank/local-cell concatenation used by the writer.
        """
        return self._ensure_cell_index()

    def _get_cellid_base(self) -> int:
        """Return the original package's one-based grid convention.

        This is intentionally independent of the values stored in the file's
        ``CellID`` array.  RHybrid output can contain zero-based CellIDs, while
        the original package's coordinate and neighbour routines define
        geometric cell IDs as 1..N.  Keeping those concepts separate is
        required to reproduce existing RHybrid analyses.
        """
        return 1

    def _invalid_cellid(self) -> int:
        return 0

    def is_valid_cellid(self, cellid: Any) -> bool | np.ndarray:
        """Return whether CellIDs are stored in the SpatialGrid array."""
        scalar = np.isscalar(cellid)
        ids = np.atleast_1d(cellid).astype(np.int64)
        index = self._ensure_cell_index()
        valid = np.array([int(value) in index for value in ids], dtype=bool)
        return bool(valid[0]) if scalar else valid

    def get_cell_indices(
        self, cellids: Any, reflevels: Any | None = None
    ) -> np.ndarray:
        """Convert uniform-grid CellIDs to zero-based ``[i, j, k]`` indices."""
        scalar = np.isscalar(cellids)
        ids = np.atleast_1d(cellids).astype(np.int64)
        if reflevels is not None:
            levels = np.broadcast_to(np.asarray(reflevels), ids.shape)
            if np.any(levels != 0):
                raise NotImplementedError(
                    "RHybridPy supports uniform refinement level 0 only"
                )
        count = int(np.prod(self._mesh_size))
        base = self._get_cellid_base()
        valid = (ids >= base) & (ids < base + count)
        zero_based = ids - base
        indices = np.full((len(ids), 3), -1, dtype=np.int64)
        nx, ny = int(self._mesh_size[0]), int(self._mesh_size[1])
        indices[valid, 0] = zero_based[valid] % nx
        indices[valid, 1] = (zero_based[valid] // nx) % ny
        indices[valid, 2] = zero_based[valid] // (nx * ny)
        return indices[0] if scalar else indices

    def get_cell_coordinates(self, cellids: Any) -> np.ndarray:
        """Return physical cell-centre coordinates for uniform-grid CellIDs."""
        scalar = np.isscalar(cellids)
        indices = np.atleast_2d(self.get_cell_indices(cellids))
        spacing = (self._extent[3:] - self._extent[:3]) / self._mesh_size
        coordinates = self._extent[:3] + (indices + 0.5) * spacing
        coordinates[np.any(indices < 0, axis=1)] = np.nan
        return coordinates[0] if scalar else coordinates

    def get_cell_dx(self, cellid: Any) -> np.ndarray:
        """Return the constant ``[dx, dy, dz]`` of one or more cells."""
        scalar = np.isscalar(cellid)
        ids = np.atleast_1d(cellid)
        spacing = (self._extent[3:] - self._extent[:3]) / self._mesh_size
        result = np.repeat(spacing[np.newaxis, :], len(ids), axis=0)
        return result[0] if scalar else result

    def get_cell_bbox(self, cellid: Any) -> tuple[np.ndarray, np.ndarray]:
        """Return lower and upper physical corners of one or more cells."""
        half = self.get_cell_dx(cellid) * 0.5
        centre = self.get_cell_coordinates(cellid)
        return centre - half, centre + half

    def get_vertex_indices(self, coordinates: Any) -> tuple[int, int, int] | list[tuple[int, int, int]]:
        """Convert physical vertex coordinates to uniform-grid vertex indices."""
        points = np.asarray(coordinates, dtype=float)
        scalar = points.ndim == 1
        points = np.atleast_2d(points)
        if points.shape[1] != 3:
            raise ValueError("coordinates must have shape (3,) or (N, 3)")
        spacing = (self._extent[3:] - self._extent[:3]) / self._mesh_size
        tolerance = np.finfo(float).eps * np.maximum(1.0, np.abs(spacing)) * 8
        indices = np.floor(
            (points - self._extent[:3] + tolerance) / spacing
        ).astype(np.int64)
        tuples = [tuple(int(value) for value in row) for row in indices]
        return tuples[0] if scalar else tuples

    def get_vertex_coordinates_from_indices(self, indices: Any) -> np.ndarray:
        """Convert uniform-grid vertex indices to physical coordinates."""
        array = np.asarray(indices, dtype=np.int64)
        scalar = array.ndim == 1
        array = np.atleast_2d(array)
        if array.shape[1] != 3:
            raise ValueError("indices must have shape (3,) or (N, 3)")
        spacing = (self._extent[3:] - self._extent[:3]) / self._mesh_size
        coordinates = self._extent[:3] + array * spacing
        return coordinates[0] if scalar else coordinates

    def get_cell_corner_vertices(
        self, cids: Any
    ) -> dict[int, tuple[tuple[int, int, int], ...]]:
        """Return the eight corner-vertex index tuples for each CellID."""
        ids = np.atleast_1d(cids).astype(np.int64)
        cell_indices = np.atleast_2d(self.get_cell_indices(ids))
        result: dict[int, tuple[tuple[int, int, int], ...]] = {}
        for cellid, low in zip(ids, cell_indices):
            if np.any(low < 0):
                raise ValueError(f"CellID {int(cellid)} is outside the uniform mesh")
            corners = tuple(
                tuple(int(value) for value in low + (dx, dy, dz))
                for dx in (0, 1)
                for dy in (0, 1)
                for dz in (0, 1)
            )
            result[int(cellid)] = corners
        return result

    def get_cell_neighbor(
        self,
        cellidss: Any,
        offsetss: Any,
        periodic: Sequence[bool],
        prune_uniques: bool = False,
    ) -> int | np.ndarray:
        """Return uniform-grid neighbors, or the file's invalid-ID sentinel."""
        del prune_uniques  # accepted for API compatibility; no cache to optimize
        scalar = np.isscalar(cellidss)
        ids = np.atleast_1d(cellidss).astype(np.int64)
        offsets = np.asarray(offsetss, dtype=np.int64)
        if offsets.ndim == 1:
            offsets = np.broadcast_to(offsets, (len(ids), 3)).copy()
        if offsets.shape != (len(ids), 3):
            raise ValueError("offsetss must have shape (3,) or (N, 3)")
        if len(periodic) != 3:
            raise ValueError("periodic must contain three booleans")

        indices = np.atleast_2d(self.get_cell_indices(ids))
        neighbors = indices + offsets
        valid = np.all(indices >= 0, axis=1)
        for axis in range(3):
            if periodic[axis]:
                neighbors[:, axis] %= self._mesh_size[axis]
            else:
                valid &= (neighbors[:, axis] >= 0) & (
                    neighbors[:, axis] < self._mesh_size[axis]
                )
        base = self._get_cellid_base()
        invalid = self._invalid_cellid()
        result = np.full(len(ids), invalid, dtype=np.int64)
        nx, ny = int(self._mesh_size[0]), int(self._mesh_size[1])
        result[valid] = (
            neighbors[valid, 0]
            + neighbors[valid, 1] * nx
            + neighbors[valid, 2] * nx * ny
            + base
        )
        available = self._ensure_cell_index()
        result = np.array(
            [cellid if int(cellid) in available else invalid for cellid in result],
            dtype=np.int64,
        )
        return int(result[0]) if scalar else result

    def get_cellid(self, coordinates: Any) -> int | np.ndarray:
        points = np.asarray(coordinates, dtype=float)
        scalar = points.ndim == 1
        points = np.atleast_2d(points)
        if points.shape[1] != 3:
            raise ValueError("coordinates must have shape (3,) or (N, 3)")
        low, high = self._extent[:3], self._extent[3:]
        # Match the original package: points exactly on a lower or upper boundary are
        # outside the mesh and receive the invalid CellID 0.
        inside = np.all((points > low) & (points < high), axis=1)
        base = self._get_cellid_base()
        invalid = self._invalid_cellid()
        result = np.full(len(points), invalid, dtype=np.int64)
        indices = np.floor(
            (points[inside] - low) / ((high - low) / self._mesh_size)
        ).astype(np.int64)
        result[inside] = (
            indices[:, 0]
            + indices[:, 1] * self._mesh_size[0]
            + indices[:, 2] * self._mesh_size[0] * self._mesh_size[1]
            + base
        )
        available = self._ensure_cell_index()
        result = np.array(
            [cid if int(cid) in available else invalid for cid in result],
            dtype=np.int64,
        )
        return int(result[0]) if scalar else result

    def read_interpolated_variable(
        self,
        name: str,
        coordinates: Any,
        operator: str | int = "pass",
        periodic: Sequence[bool] = (True, True, True),
        method: str = "linear",
    ) -> Any:
        """Interpolate a SpatialGrid variable on a uniform cell-centred mesh."""
        if method.lower() not in {"linear", "nearest"}:
            raise ValueError("method must be 'linear' or 'nearest'")
        points = np.asarray(coordinates, dtype=float)
        scalar = points.ndim == 1
        points = np.atleast_2d(points)
        if points.shape[1] != 3:
            raise ValueError("coordinates must have shape (3,) or (N, 3)")

        closest = np.asarray(self.get_cellid(points), dtype=np.int64)
        sample = np.asarray(self.read_variable(name, operator=operator))
        value_shape = sample.shape[1:] if sample.ndim > 1 else ()

        if method.lower() == "nearest":
            result = np.full(
                (len(points),) + value_shape,
                np.nan,
                dtype=np.result_type(sample, float),
            )
            valid = closest != 0
            if np.any(valid):
                result[valid] = self.read_variable(
                    name, cellids=closest[valid], operator=operator
                )
            return result[0] if scalar else result

        result = np.full(
            (len(points),) + value_shape,
            np.nan,
            dtype=np.result_type(sample, float),
        )
        valid_points = closest != 0
        if not np.any(valid_points):
            return result[0] if scalar else result

        use_points = points[valid_points]
        use_closest = closest[valid_points]
        closest_coordinates = np.atleast_2d(self.get_cell_coordinates(use_closest))
        lower_offsets = np.zeros((len(use_points), 3), dtype=np.int64)
        lower_offsets[use_points <= closest_coordinates] = -1
        lower_ids = np.atleast_1d(
            self.get_cell_neighbor(use_closest, lower_offsets, periodic)
        ).astype(np.int64)

        corner_offsets = np.array(
            [
                (x, y, z)
                for x in (0, 1)
                for y in (0, 1)
                for z in (0, 1)
            ],
            dtype=np.int64,
        )
        neighbor_ids = self.get_cell_neighbor(
            np.repeat(lower_ids, 8),
            np.tile(corner_offsets, (len(lower_ids), 1)),
            periodic,
        ).reshape(len(lower_ids), 8)

        neighbor_values = np.full(
            (len(lower_ids), 8) + value_shape,
            np.nan,
            dtype=np.result_type(sample, float),
        )
        flat_ids = neighbor_ids.reshape(-1)
        stored = self._ensure_cell_index()
        readable = np.array([int(cid) in stored for cid in flat_ids])
        if np.any(readable):
            flat_values = neighbor_values.reshape((-1,) + value_shape)
            flat_values[readable] = self.read_variable(
                name, cellids=flat_ids[readable], operator=operator
            )
        neighbor_values = neighbor_values.reshape(
            (len(lower_ids), 2, 2, 2) + value_shape
        )

        lower_coordinates = np.atleast_2d(self.get_cell_coordinates(lower_ids))
        upper_coordinates = np.atleast_2d(
            self.get_cell_coordinates(neighbor_ids[:, 7])
        )
        scaled = np.zeros_like(use_points)
        changing = lower_coordinates != upper_coordinates
        scaled[changing] = (
            (use_points - lower_coordinates)[changing]
            / (upper_coordinates - lower_coordinates)[changing]
        )

        extra = tuple(range(3, 3 + len(value_shape)))
        c2 = neighbor_values[:, 0] * (
            1 - np.expand_dims(scaled[:, 0], axis=(1, 2) + extra)
        ) + neighbor_values[:, 1] * np.expand_dims(
            scaled[:, 0], axis=(1, 2) + extra
        )
        extra = tuple(range(2, 2 + len(value_shape)))
        c1 = c2[:, 0] * (
            1 - np.expand_dims(scaled[:, 1], axis=(1,) + extra)
        ) + c2[:, 1] * np.expand_dims(scaled[:, 1], axis=(1,) + extra)
        extra = tuple(range(1, 1 + len(value_shape)))
        interpolated = c1[:, 0] * (
            1 - np.expand_dims(scaled[:, 2], axis=extra)
        ) + c1[:, 1] * np.expand_dims(scaled[:, 2], axis=extra)
        result[valid_points] = interpolated
        return result[0] if scalar else result


def _is_all_cells(cellids: Any) -> bool:
    return np.isscalar(cellids) and int(cellids) == -1


class VlsvParticles:
    """Read point meshes and point-data variables from a VLSV file."""

    def __init__(self, file_name: str) -> None:
        self.file_name = os.path.abspath(os.fspath(file_name))
        if not os.path.isfile(self.file_name):
            raise FileNotFoundError(self.file_name)
        self.__xml_root = _read_xml_footer(self.file_name)

    def _read(self, child: ET.Element) -> np.ndarray:
        vector_size = int(_literal(child.attrib.get("vectorsize"), 1))
        array_size = int(_literal(child.attrib.get("arraysize"), 1))
        element_size = int(_literal(child.attrib.get("datasize")))
        offset = int(_literal(child.text))
        with open(self.file_name, "rb") as stream:
            stream.seek(offset)
            data = np.fromfile(
                stream,
                dtype=_dtype(child.attrib["datatype"], element_size),
                count=array_size * vector_size,
            )
        if data.size != array_size * vector_size:
            raise IOError(f"Short particle-data read for {child.attrib.get('name')!r}")
        return data.reshape(array_size, vector_size)

    def get_all_point_meshes(self) -> list[str]:
        return [
            child.attrib["name"]
            for child in self.__xml_root
            if child.tag == "MESH"
            and child.attrib.get("type") == "point"
            and "name" in child.attrib
        ]

    def get_all_point_variables(self) -> list[str]:
        return [
            child.attrib["name"]
            for child in self.__xml_root
            if child.tag == "VARIABLE"
            and child.attrib.get("type") == "pointdata"
            and "name" in child.attrib
        ]

    def read_particle_mesh(self, name: str) -> np.ndarray | None:
        for child in self.__xml_root:
            if (
                child.tag == "MESH"
                and child.attrib.get("type") == "point"
                and child.attrib.get("name") == name
            ):
                return self._read(child)
        return None

    def read_particle_variable(self, name: str) -> np.ndarray | None:
        for child in self.__xml_root:
            if (
                child.tag == "VARIABLE"
                and child.attrib.get("type") == "pointdata"
                and child.attrib.get("name") == name
            ):
                return self._read(child)
        return None


class VlsvWriter:
    """Write basic VLSV variables using an existing reader as the template.

    New bulk arrays must use the template reader's file order.  Arrays derived
    elementwise from variables read from that same reader already have the
    correct order because the copied ``CellID`` array preserves the mapping.
    """

    _COPY_TAGS = {
        "PARAMETER",
        "PARAMETERS",
        "MESH_NODE_CRDS",
        "MESH_NODE_CRDS_X",
        "MESH_NODE_CRDS_Y",
        "MESH_NODE_CRDS_Z",
        "MESH_OFFSETS",
        "MESH",
        "MESH_DOMAIN_SIZES",
        "MESH_DOMAIN_EXTENTS",
        "MESH_GHOST_DOMAINS",
        "MESH_GHOST_LOCALIDS",
        "MESH_BBOX",
        "COORDS",
    }

    def __init__(
        self,
        vlsvReader: VlsvReader,
        file_name: str,
        copy_meshes: Sequence[str] | None = None,
        clone: bool = False,
    ) -> None:
        self.file_name = os.path.abspath(os.fspath(file_name))
        self.__xml_root = ET.Element("VLSV")
        self.__closed = False
        self._written_entries: set[tuple[str, str, str]] = set()
        if clone:
            shutil.copy2(vlsvReader.file_name, self.file_name)
            self.__xml_root = ET.fromstring(
                ET.tostring(vlsvReader._VlsvReader__xml_root)
            )
            self.__fptr = open(self.file_name, "r+b")
            self.__fptr.seek(8)
            (xml_offset,) = struct.unpack("=Q", self.__fptr.read(8))
            self.__fptr.seek(xml_offset)
            self._written_entries = {
                (
                    child.tag,
                    child.attrib.get("name", "").lower(),
                    child.attrib.get("mesh", ""),
                )
                for child in self.__xml_root
            }
            return

        self.__fptr = open(self.file_name, "w+b")
        np.array(0, dtype=np.uint64).tofile(self.__fptr)
        np.array(0, dtype=np.uint64).tofile(self.__fptr)
        self._initialize(vlsvReader, copy_meshes)

    def __enter__(self) -> "VlsvWriter":
        return self

    def __exit__(self, *_: Any) -> None:
        self.close()

    def _initialize(
        self, reader: VlsvReader, copy_meshes: Sequence[str] | None
    ) -> None:
        selected = None if copy_meshes is None else set(copy_meshes)
        selected_cellid_meshes = {
            child.attrib.get("mesh", "")
            for child in reader._VlsvReader__xml_root
            if child.tag == "VARIABLE"
            and child.attrib.get("name", "").lower() == "cellid"
            and child.attrib.get("mesh", "")
            and (selected is None or child.attrib.get("mesh") in selected)
        }
        for child in reader._VlsvReader__xml_root:
            is_cellid = (
                child.tag == "VARIABLE"
                and child.attrib.get("name", "").lower() == "cellid"
            )
            if child.tag not in self._COPY_TAGS and not is_cellid:
                continue
            mesh = child.attrib.get("mesh", "")
            name = child.attrib.get("name", "")
            if selected is not None:
                if mesh and mesh not in selected:
                    continue
                if child.tag == "MESH" and name not in selected:
                    continue
            if is_cellid and not mesh and selected_cellid_meshes:
                continue
            extra = {
                key: value
                for key, value in child.attrib.items()
                if key not in {"name", "mesh", "arraysize", "vectorsize", "datatype", "datasize"}
            }
            data = reader._read_element(child)
            self._write(data, name=name, tag=child.tag, mesh=mesh, extra_attribs=extra)

    def copy_variables(
        self, vlsvReader: VlsvReader, varlist: Sequence[str] | None = None
    ) -> None:
        wanted = None if varlist is None else {name.lower() for name in varlist}
        found: set[str] = set()
        for child in vlsvReader._VlsvReader__xml_root:
            if child.tag != "VARIABLE" or "name" not in child.attrib:
                continue
            name = child.attrib["name"]
            if wanted is not None and name.lower() not in wanted:
                continue
            mesh = child.attrib.get("mesh", "")
            entry = ("VARIABLE", name.lower(), mesh)
            if entry in self._written_entries:
                found.add(name.lower())
                continue
            extra = {
                key: value
                for key, value in child.attrib.items()
                if key not in {"name", "mesh", "arraysize", "vectorsize", "datatype", "datasize"}
            }
            self._write(
                vlsvReader._read_element(child),
                name=name,
                tag="VARIABLE",
                mesh=mesh,
                extra_attribs=extra,
            )
            found.add(name.lower())
        if wanted is not None:
            missing = wanted - found
            if missing:
                raise KeyError(
                    "Variables not stored in the source file: "
                    + ", ".join(sorted(missing))
                )

    def write(
        self,
        data: Any,
        name: str,
        tag: str,
        mesh: str,
        extra_attribs: dict[str, Any] | None = None,
    ) -> bool:
        if tag == "VARIABLE":
            warnings.warn("Prefer write_variable_info() for variables", stacklevel=2)
        return self._write(data, name, tag, mesh, extra_attribs or {})

    def _write(
        self,
        data: Any,
        name: str,
        tag: str,
        mesh: str | None,
        extra_attribs: dict[str, Any] | None = None,
    ) -> bool:
        if self.__closed:
            raise ValueError("I/O operation on closed VlsvWriter")
        if data is None:
            warnings.warn(f"Skipping None data for {name!r}", stacklevel=2)
            return False
        array = np.atleast_1d(np.ma.getdata(data))
        if array.ndim > 2:
            raise ValueError("VLSV arrays must be one- or two-dimensional")
        if array.dtype.kind not in "fiu" or array.dtype.itemsize not in (4, 8):
            if array.dtype.kind in "iu":
                array = array.astype(np.int64 if array.dtype.kind == "i" else np.uint64)
            elif array.dtype.kind == "f":
                array = array.astype(np.float64)
            else:
                raise TypeError(f"Unsupported array dtype {array.dtype}")
        # VLSV files are little-endian in the supported workflow.
        array = np.ascontiguousarray(array.astype(array.dtype.newbyteorder("<"), copy=False))

        child = ET.SubElement(self.__xml_root, tag)
        child.attrib["name"] = str(name)
        if mesh:
            child.attrib["mesh"] = str(mesh)
        child.attrib["arraysize"] = str(array.shape[0])
        child.attrib["vectorsize"] = str(array.shape[1] if array.ndim == 2 else 1)
        child.attrib["datatype"] = {"f": "float", "i": "int", "u": "uint"}[array.dtype.kind]
        child.attrib["datasize"] = str(array.dtype.itemsize)
        for key, value in (extra_attribs or {}).items():
            child.attrib[str(key)] = str(value)
        child.text = str(self.__fptr.tell())
        array.tofile(self.__fptr)
        self._written_entries.add((tag, str(name).lower(), str(mesh or "")))
        self._write_xml_footer()
        return True

    def write_variable_info(
        self,
        varinfo: VariableInfo,
        mesh: str,
        unitConversion: str | int | float,
        extra_attribs: dict[str, Any] | None = None,
    ) -> bool:
        attributes = dict(extra_attribs or {})
        attributes.update(
            {
                "variableLaTeX": varinfo.latex,
                "unit": varinfo.units,
                "unitLaTeX": varinfo.latexunits,
                "unitConversion": unitConversion,
            }
        )
        return self._write(
            varinfo.data,
            name=varinfo.name,
            tag="VARIABLE",
            mesh=mesh,
            extra_attribs=attributes,
        )

    def write_variable(
        self,
        data: Any,
        name: str,
        mesh: str,
        units: str = "",
        latex: str = "",
        latexunits: str = "",
        unitConversion: str | int | float = 1,
        extra_attribs: dict[str, Any] | None = None,
    ) -> bool:
        return self.write_variable_info(
            VariableInfo(data, name=name, units=units, latex=latex, latexunits=latexunits),
            mesh,
            unitConversion,
            extra_attribs,
        )

    def _write_xml_footer(self) -> None:
        data_end = self.__fptr.tell()
        root_copy = ET.fromstring(ET.tostring(self.__xml_root))
        _indent_xml(root_copy)
        self.__fptr.seek(data_end)
        self.__fptr.truncate()
        ET.ElementTree(root_copy).write(
            self.__fptr, encoding="utf-8", xml_declaration=False
        )
        self.__fptr.seek(8)
        np.array(data_end, dtype=np.uint64).tofile(self.__fptr)
        self.__fptr.seek(data_end)
        self.__fptr.flush()

    def close(self) -> None:
        if not self.__closed:
            self._write_xml_footer()
            self.__fptr.close()
            self.__closed = True


def _indent_xml(element: ET.Element, level: int = 0) -> None:
    indent = "\n" + level * "   "
    if len(element):
        if not element.text or not element.text.strip():
            element.text = indent + "   "
        for child in element:
            _indent_xml(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = indent
    elif level and (not element.tail or not element.tail.strip()):
        element.tail = indent


def vlsv_intpol_points(
    vlsvReader: VlsvReader,
    points: Any,
    varlist: Sequence[str],
    operator: str | int = "pass",
    interpolation_order: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Read one or more variables at points on a uniform spatial mesh."""
    coordinates = np.asarray(points, dtype=float)
    if coordinates.ndim == 1:
        coordinates = np.atleast_2d(coordinates)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("points must have shape (3,) or (N, 3)")
    if interpolation_order not in (0, 1):
        raise ValueError("interpolation_order must be 0 or 1")
    names = list(varlist) or vlsvReader.get_all_variables()
    if not names:
        raise ValueError("No variables requested or stored")

    cellids = np.asarray(vlsvReader.get_cellid(coordinates)).reshape(-1, 1)
    columns: list[np.ndarray] = []
    header = "x y z cellid "
    for name in names:
        if not vlsvReader.check_variable(name):
            raise KeyError(
                f"Variable {name!r} does not exist in {vlsvReader.file_name!r}"
            )
        if interpolation_order == 1:
            values = vlsvReader.read_interpolated_variable(
                name, coordinates, operator=operator
            )
        else:
            sample = np.asarray(vlsvReader.read_variable(name, operator=operator))
            value_shape = sample.shape[1:] if sample.ndim > 1 else ()
            values = np.full(
                (len(coordinates),) + value_shape, np.nan, dtype=np.result_type(sample, float)
            )
            # The original package reserves geometric CellID 0 as the out-of-domain
            # sentinel even when an RHybrid file stores a data key named 0.
            valid = cellids[:, 0] != 0
            if np.any(valid):
                values[valid] = vlsvReader.read_variable(
                    name, cellids=cellids[valid, 0], operator=operator
                )
        values = np.asarray(values)
        if values.ndim == 1:
            values = values[:, None]
        columns.append(values.reshape(len(coordinates), -1))
        if values.shape[1] == 1:
            header += f"{name} "
        else:
            header += "".join(f"{name}({i}) " for i in range(values.shape[1]))
    return coordinates, cellids, np.concatenate(columns, axis=1), header


# Namespace shims retain the familiar ``rhb.vlsvfile`` and
# ``rhb.calculations`` API while keeping the implementation in this one file.
vlsvfile = types.SimpleNamespace(
    VlsvReader=VlsvReader,
    VlsvParticles=VlsvParticles,
    VlsvWriter=VlsvWriter,
)
calculations = types.SimpleNamespace(
    VariableInfo=VariableInfo,
    vlsv_intpol_points=vlsv_intpol_points,
)
