#!/usr/bin/env python3
"""Average SpatialGrid bulk variables from multiple VLSV files.

The first input file is the output template, and its CellID storage order
defines the output order.  Every input is aligned by CellID before arithmetic,
so MPI load-balancing changes between timesteps do not mix values belonging to
different cells.  The streaming strategy is intended for large grids: it uses
balanced variable groups, transient uint32 CellID maps, chunked reordering,
and temporary float64 memory-mapped accumulators.

Examples
--------
    python average_state_files.py -o average.vlsv state*.vlsv

    python average_state_files.py -o average.vlsv \
        --variables cellB,ne,cellUe state0000*.vlsv

    python average_state_files.py -o average.vlsv \
        --strategy file-major state*.vlsv

    python average_state_files.py -o average.vlsv \
        --strategy variable-major --workers 8 state*.vlsv

    python average_state_files.py -o average.vlsv \
        --strategy streaming --workers 6 --work-dir "$SLURM_TMPDIR" state*.vlsv

    python average_state_files.py -o average.vlsv \
        --strategy streaming --workers 6 --input-list input_files.txt

The output inherits mesh information and parameters (including the time
parameter) from the first input.  Since a temporal average has no unique time,
interpret or adjust that copied time metadata according to your workflow.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import getpass
import hashlib
import multiprocessing as mp
import os
from pathlib import Path
import shlex
import shutil
import socket
import sys
import tempfile
import time
from typing import Iterable, Sequence

import numpy as np

import rhybridpy as rhb


_STRUCTURAL_ATTRIBUTES = {
    "name",
    "mesh",
    "arraysize",
    "vectorsize",
    "datatype",
    "datasize",
    "unit",
    "unitLaTeX",
    "variableLaTeX",
    "unitConversion",
}


# Initialized separately inside every multiprocessing worker.  Workers only
# read VLSV files; the parent process remains the sole output writer.
_WORKER_READERS: list[rhb.VlsvReader] | None = None
_WORKER_ALIGNMENTS: Sequence[np.ndarray | None] | None = None
_WORKER_PATHS: Sequence[str] | None = None

# State used only by the low-memory streaming worker strategy.
_STREAM_PATHS: Sequence[str] | None = None
_STREAM_REFERENCE_KEYS: np.ndarray | None = None
_STREAM_MINIMUM_CELLID: int | None = None
_STREAM_CELL_COUNT: int | None = None
_STREAM_CHUNK_CELLS: int | None = None
_STREAM_RAW_POSITIONS: np.ndarray | None = None
_STREAM_POSITIONS: np.ndarray | None = None


def _progress_step(total: int, updates: int = 10) -> int:
    """Return an interval producing at most approximately ``updates`` lines."""
    return max(1, (total + updates - 1) // updates)


def _report_progress(index: int, total: int, step: int) -> bool:
    return index == total or index == 1 or index % step == 0


def _format_gib(byte_count: int) -> str:
    return f"{byte_count / 1024**3:.3f} GiB"


def _markdown_code(value: object) -> str:
    """Return a value protected for a Markdown inline-code span."""
    return str(value).replace("`", r"\`")


def _format_launch_command() -> list[str]:
    """Return the exact Python invocation with an optional launcher prefix."""
    prefix = os.environ.get("AVERAGE_STATE_LAUNCH_PREFIX")
    if prefix:
        tokens = [*shlex.split(prefix), *sys.argv]
    else:
        tokens = [sys.executable, *sys.argv]

    lines = ["```bash"]
    for index, token in enumerate(tokens):
        continuation = " \\" if index + 1 < len(tokens) else ""
        indentation = "" if index == 0 else "    "
        lines.append(f"{indentation}{shlex.quote(token)}{continuation}")
    lines.append("```")
    return lines


def _suggested_snapshot_context(
    input_paths: Sequence[Path],
) -> tuple[list[str], list[Path]]:
    """Return Markdown selection lines and absent suggested snapshot paths."""
    names = {
        "first": "AVERAGE_STATE_TIMESTEP_FIRST",
        "last": "AVERAGE_STATE_TIMESTEP_LAST",
        "stride": "AVERAGE_STATE_TIMESTEP_STRIDE",
        "width": "AVERAGE_STATE_NUMBER_WIDTH",
        "folder": "AVERAGE_STATE_INPUT_FOLDER",
    }
    values = {name: os.environ.get(variable) for name, variable in names.items()}
    if any(value is None for value in values.values()):
        return [f"- Input files supplied: {len(input_paths):,}"], []

    try:
        first = int(values["first"] or "")
        last = int(values["last"] or "")
        stride = int(values["stride"] or "")
        width = int(values["width"] or "")
    except ValueError:
        return [
            f"- Input files supplied: {len(input_paths):,}",
            "- Suggested-range metadata: invalid and therefore omitted",
        ], []

    if first < 0 or last < first or stride < 1 or width < 1:
        return [
            f"- Input files supplied: {len(input_paths):,}",
            "- Suggested-range metadata: invalid and therefore omitted",
        ], []

    folder = Path(values["folder"] or ".").expanduser().resolve()
    suggested = [
        (folder / f"state{timestep:0{width}d}.vlsv").resolve()
        for timestep in range(first, last + 1, stride)
    ]
    actual = set(input_paths)
    missing = [path for path in suggested if path not in actual]
    lines = [
        f"- Suggested first timestep: {first}",
        f"- Suggested last timestep: {last}",
        f"- Suggested stride: {stride}",
        f"- State-number width: {width}",
        f"- Suggested files: {len(suggested):,}",
        f"- Found files: {len(input_paths):,}",
        f"- Missing files: {len(missing):,}",
    ]
    return lines, missing


def _write_run_log(
    *,
    output: Path,
    input_paths: Sequence[Path],
    selected: Sequence[str],
    cell_count: int,
    strategy: str,
    requested_workers: int,
    effective_workers: int,
    chunk_cells: int,
    started_utc: datetime,
    elapsed_seconds: float,
    overwrite: bool,
) -> Path:
    """Write a human-readable Markdown provenance log beside the output."""
    log_path = Path(f"{output}.log")
    if log_path.exists() and not overwrite:
        raise FileExistsError(f"Run log already exists: {log_path}; use --force")

    finished_utc = datetime.now(timezone.utc)
    selection_lines, missing_paths = _suggested_snapshot_context(input_paths)
    script_path = Path(__file__).resolve()
    script_hash = hashlib.sha256(script_path.read_bytes()).hexdigest()
    loaded_modules = os.environ.get("LOADEDMODULES") or "not reported"
    rhybridpy_version = getattr(rhb, "__version__", "not reported")

    lines = [
        "# VLSV state-file average",
        "",
        "## Result",
        "",
        f"- Output file: `{_markdown_code(output)}`",
        f"- Output size: {output.stat().st_size:,} bytes",
        f"- Started: `{started_utc.isoformat(timespec='seconds')}`",
        f"- Finished: `{finished_utc.isoformat(timespec='seconds')}`",
        f"- Elapsed: {elapsed_seconds:.1f} seconds",
        f"- User: `{_markdown_code(getpass.getuser())}`",
        f"- Host: `{_markdown_code(socket.gethostname())}`",
        f"- Working directory: `{_markdown_code(Path.cwd())}`",
        f"- Slurm job ID: `{_markdown_code(os.environ.get('SLURM_JOB_ID', 'not applicable'))}`",
        f"- Slurm cluster: `{_markdown_code(os.environ.get('SLURM_CLUSTER_NAME', 'not reported'))}`",
        f"- Slurm partition: `{_markdown_code(os.environ.get('SLURM_JOB_PARTITION', 'not reported'))}`",
        f"- Slurm nodes: `{_markdown_code(os.environ.get('SLURM_NNODES', 'not reported'))}`",
        f"- Slurm tasks: `{_markdown_code(os.environ.get('SLURM_NTASKS', 'not reported'))}`",
        f"- Slurm CPUs per task: `{_markdown_code(os.environ.get('SLURM_CPUS_PER_TASK', 'not reported'))}`",
        f"- Slurm memory per node: `{_markdown_code(os.environ.get('SLURM_MEM_PER_NODE', 'not reported'))}`",
        "",
        "## Averaging configuration",
        "",
        f"- Strategy: `{strategy}`",
        f"- Requested workers: {requested_workers}",
        f"- Effective workers: {effective_workers}",
        f"- Chunk size: {chunk_cells:,} cells",
        f"- Grid cells: {cell_count:,}",
        f"- Selected variables: {len(selected):,}",
        "",
        "### Variables",
        "",
        *[f"- `{_markdown_code(name)}`" for name in selected],
        "",
        "## Snapshot selection",
        "",
        *selection_lines,
        "",
        "## Software environment",
        "",
        f"- Averaging script: `{_markdown_code(script_path)}`",
        f"- Averaging-script SHA-256: `{script_hash}`",
        f"- Python: `{_markdown_code(sys.version.split()[0])}`",
        f"- NumPy: `{_markdown_code(np.__version__)}`",
        f"- RHybridPy: `{_markdown_code(rhybridpy_version)}`",
        f"- Loaded modules: `{_markdown_code(loaded_modules)}`",
        "",
        "## Launch command",
        "",
        *_format_launch_command(),
    ]

    if missing_paths:
        lines.extend(
            [
                "",
                "## Missing suggested files",
                "",
                *[
                    f"{index}. `{_markdown_code(path)}`"
                    for index, path in enumerate(missing_paths, start=1)
                ],
            ]
        )

    # Keep the potentially long actual input list at the very end.
    lines.extend(
        [
            "",
            "## Input files used",
            "",
            *[
                f"{index}. `{_markdown_code(path)}`"
                for index, path in enumerate(input_paths, start=1)
            ],
            "",
        ]
    )

    temporary_log = log_path.with_name(f".{log_path.name}.{os.getpid()}.tmp")
    try:
        temporary_log.write_text("\n".join(lines), encoding="utf-8")
        os.replace(temporary_log, log_path)
    finally:
        temporary_log.unlink(missing_ok=True)
    return log_path


def _bulk_variable_elements(reader: rhb.VlsvReader) -> dict[str, object]:
    """Return one XML element per readable SpatialGrid bulk variable."""
    exact: dict[str, object] = {}
    legacy: dict[str, object] = {}

    for child in reader.xml_root:
        if child.tag != "VARIABLE" or "name" not in child.attrib:
            continue
        name = child.attrib["name"]
        if name.lower() == "cellid" or child.attrib.get("type") == "pointdata":
            continue
        mesh = child.attrib.get("mesh")
        if mesh == "SpatialGrid":
            exact.setdefault(name, child)
        elif mesh is None:
            # RHybridPy supports old meshless bulk entries as a fallback.
            legacy.setdefault(name, child)

    # Exact SpatialGrid entries take precedence over legacy meshless entries.
    return {**legacy, **exact}


def _parse_variables(value: str | None) -> list[str] | None:
    if value is None:
        return None
    names = [name.strip() for name in value.split(",") if name.strip()]
    if not names:
        raise argparse.ArgumentTypeError("--variables must contain a name")
    if len(names) != len(set(names)):
        raise argparse.ArgumentTypeError("--variables contains duplicates")
    return names


def _read_input_list(path: os.PathLike[str] | str) -> list[str]:
    """Read one input path per line without imposing a command-size limit."""
    list_path = Path(path).expanduser().resolve()
    if not list_path.is_file():
        raise FileNotFoundError(f"Input-list file not found: {list_path}")
    paths = [
        line
        for line in list_path.read_text(encoding="utf-8").splitlines()
        if line
    ]
    if not paths:
        raise ValueError(f"Input-list file is empty: {list_path}")
    return paths


def _validate_inputs(
    readers: Sequence[rhb.VlsvReader],
    input_paths: Sequence[Path],
    variables: Sequence[str] | None,
) -> tuple[list[np.ndarray | None], list[str], dict[str, object]]:
    reference = readers[0]
    reference_elements = _bulk_variable_elements(reference)

    if variables is None:
        selected = sorted(reference_elements)
    else:
        missing = [name for name in variables if name not in reference_elements]
        if missing:
            raise ValueError(
                f"Variables absent from first input: {', '.join(missing)}"
            )
        selected = list(variables)

    if not selected:
        raise ValueError("No SpatialGrid bulk variables were selected")
    reference_size = reference.get_spatial_mesh_size()
    reference_extent = reference.get_spatial_mesh_extent()
    reference_cellids = np.asarray(reference.read_variable("CellID"), dtype=np.int64)
    if reference_cellids.ndim != 1:
        raise ValueError("CellID must be a one-dimensional array")
    reference_sort_order = np.argsort(reference_cellids, kind="stable")
    reference_id_set = reference_cellids[reference_sort_order]
    if np.any(reference_id_set[1:] == reference_id_set[:-1]):
        raise ValueError(f"Duplicate CellIDs in {input_paths[0]}")

    selected_set = set(selected)
    # None means that a file already has the reference/output storage order.
    alignments: list[np.ndarray | None] = [None]

    validation_start = time.perf_counter()
    validation_step = _progress_step(len(readers))
    maximum_alignment_bytes = (
        max(0, len(readers) - 1) * reference_cellids.size * np.dtype(np.int64).itemsize
    )
    print(
        "Phase 2/3: validating grids, variables, and CellIDs; "
        f"worst-case alignment storage {_format_gib(maximum_alignment_bytes)}",
        flush=True,
    )
    print(f"  validated 1/{len(readers)}: {input_paths[0].name}", flush=True)

    for file_index, (reader, path) in enumerate(
        zip(readers[1:], input_paths[1:]), start=2
    ):
        if not np.array_equal(reader.get_spatial_mesh_size(), reference_size):
            raise ValueError(f"SpatialGrid size differs in {path}")
        if not np.array_equal(reader.get_spatial_mesh_extent(), reference_extent):
            raise ValueError(f"SpatialGrid extent differs in {path}")

        elements = _bulk_variable_elements(reader)
        if variables is None and set(elements) != set(reference_elements):
            absent = sorted(set(reference_elements) - set(elements))
            extra = sorted(set(elements) - set(reference_elements))
            raise ValueError(
                f"Bulk-variable set differs in {path}; "
                f"missing={absent}, extra={extra}"
            )
        missing = sorted(selected_set - set(elements))
        if missing:
            raise ValueError(f"Variables missing from {path}: {', '.join(missing)}")

        cellids = np.asarray(reader.read_variable("CellID"), dtype=np.int64)
        if cellids.ndim != 1 or cellids.size != reference_cellids.size:
            raise ValueError(f"CellID array size differs in {path}")
        sort_order = np.argsort(cellids, kind="stable")
        sorted_cellids = cellids[sort_order]
        if np.any(sorted_cellids[1:] == sorted_cellids[:-1]):
            raise ValueError(f"Duplicate CellIDs in {path}")
        if not np.array_equal(sorted_cellids, reference_id_set):
            raise ValueError(f"CellID set differs in {path}")

        if np.array_equal(cellids, reference_cellids):
            alignments.append(None)
        else:
            # searchsorted locates each reference CellID in the sorted input;
            # sort_order translates those locations to raw input-array indices.
            alignments.append(
                sort_order[np.searchsorted(sorted_cellids, reference_cellids)]
            )

        if _report_progress(file_index, len(readers), validation_step):
            elapsed = time.perf_counter() - validation_start
            print(
                f"  validated {file_index}/{len(readers)}: {path.name} "
                f"({elapsed:.1f} s)",
                flush=True,
            )

    alignment_bytes = sum(
        alignment.nbytes for alignment in alignments if alignment is not None
    )
    reordered_files = sum(alignment is not None for alignment in alignments)
    print(
        f"Phase 2/3 complete: {reordered_files}/{len(readers)} files require "
        f"reordering; stored alignments {_format_gib(alignment_bytes)}; "
        f"elapsed {time.perf_counter() - validation_start:.1f} s",
        flush=True,
    )

    return alignments, selected, reference_elements


def _validate_streaming_inputs(
    readers: Sequence[rhb.VlsvReader],
    input_paths: Sequence[Path],
    variables: Sequence[str] | None,
) -> tuple[np.ndarray, int, list[str], dict[str, object]]:
    """Validate metadata without retaining a CellID permutation per file.

    Streaming uses a linear-time dense CellID lookup.  The stored CellIDs may
    start at any non-negative value, but they must form one contiguous range.
    """
    reference = readers[0]
    reference_elements = _bulk_variable_elements(reference)
    if variables is None:
        selected = sorted(reference_elements)
    else:
        missing = [name for name in variables if name not in reference_elements]
        if missing:
            raise ValueError(
                f"Variables absent from first input: {', '.join(missing)}"
            )
        selected = list(variables)
    if not selected:
        raise ValueError("No SpatialGrid bulk variables were selected")
    legacy = [
        name
        for name in selected
        if reference_elements[name].attrib.get("mesh") != "SpatialGrid"
    ]
    if legacy:
        raise ValueError(
            "Streaming strategy requires explicitly SpatialGrid-tagged variables; "
            "use variable-major for legacy meshless variables: "
            + ", ".join(legacy)
        )

    reference_size = reference.get_spatial_mesh_size()
    reference_extent = reference.get_spatial_mesh_extent()
    selected_set = set(selected)
    validation_start = time.perf_counter()
    validation_step = _progress_step(len(readers))
    print(
        "Phase 2/3: validating metadata; CellID sets will be validated "
        "on demand inside streaming workers",
        flush=True,
    )
    print(f"  validated metadata 1/{len(readers)}: {input_paths[0].name}", flush=True)

    for file_index, (reader, path) in enumerate(
        zip(readers[1:], input_paths[1:]), start=2
    ):
        if not np.array_equal(reader.get_spatial_mesh_size(), reference_size):
            raise ValueError(f"SpatialGrid size differs in {path}")
        if not np.array_equal(reader.get_spatial_mesh_extent(), reference_extent):
            raise ValueError(f"SpatialGrid extent differs in {path}")

        elements = _bulk_variable_elements(reader)
        if variables is None and set(elements) != set(reference_elements):
            absent = sorted(set(reference_elements) - set(elements))
            extra = sorted(set(elements) - set(reference_elements))
            raise ValueError(
                f"Bulk-variable set differs in {path}; "
                f"missing={absent}, extra={extra}"
            )
        missing = sorted(selected_set - set(elements))
        if missing:
            raise ValueError(f"Variables missing from {path}: {', '.join(missing)}")

        if _report_progress(file_index, len(readers), validation_step):
            print(
                f"  validated metadata {file_index}/{len(readers)}: {path.name} "
                f"({time.perf_counter() - validation_start:.1f} s)",
                flush=True,
            )

    reference_cellids = np.asarray(
        reference.read_variable("CellID"), dtype=np.int64
    )
    if reference_cellids.ndim != 1 or reference_cellids.size == 0:
        raise ValueError("CellID must be a non-empty one-dimensional array")
    cell_count = int(reference_cellids.size)
    if cell_count >= int(np.iinfo(np.uint32).max):
        raise ValueError("Streaming strategy requires fewer than 2**32-1 cells")

    minimum_cellid = int(reference_cellids.min())
    maximum_cellid = int(reference_cellids.max())
    if minimum_cellid < 0 or maximum_cellid - minimum_cellid + 1 != cell_count:
        raise ValueError(
            "Streaming strategy requires dense, contiguous, non-negative "
            "stored CellIDs; use variable-major for a sparse CellID set"
        )

    reference_cellids -= minimum_cellid
    reference_keys = reference_cellids.astype(np.uint32, copy=False)
    seen = np.zeros(cell_count, dtype=np.bool_)
    seen[reference_keys] = True
    if not bool(seen.all()):
        raise ValueError(f"Duplicate or missing CellIDs in {input_paths[0]}")
    del seen, reference_cellids

    print(
        "Phase 2/3 complete: metadata valid; reference CellIDs are dense "
        f"[{minimum_cellid}, {maximum_cellid}]; elapsed "
        f"{time.perf_counter() - validation_start:.1f} s",
        flush=True,
    )
    return reference_keys, minimum_cellid, selected, reference_elements


def _read_aligned(
    reader: rhb.VlsvReader,
    name: str,
    alignment: np.ndarray | None,
) -> np.ndarray:
    """Read one variable and put it in the reference file's CellID order."""
    values = np.asarray(reader.read_variable(name))
    return values if alignment is None else values[alignment]


def _write_average(
    writer: rhb.VlsvWriter,
    elements: dict[str, object],
    name: str,
    data: np.ndarray,
) -> None:
    """Write averaged data with metadata copied from the reference variable."""
    child = elements[name]
    latex = child.attrib.get("variableLaTeX", "") or name.replace("_", r"\_")
    info = rhb.VariableInfo(
        data,
        name=name,
        units=child.attrib.get("unit", ""),
        latex=latex,
        latexunits=child.attrib.get("unitLaTeX", ""),
    )
    mesh = child.attrib.get("mesh", "")
    unit_conversion = child.attrib.get("unitConversion", 1)
    extra_attributes = {
        key: value
        for key, value in child.attrib.items()
        if key not in _STRUCTURAL_ATTRIBUTES
    }
    writer.write_variable_info(
        info,
        mesh,
        unit_conversion,
        extra_attribs=extra_attributes,
    )


def _average_file_major(
    readers: Sequence[rhb.VlsvReader],
    paths: Sequence[Path],
    alignments: Sequence[np.ndarray | None],
    selected: Sequence[str],
) -> dict[str, np.ndarray]:
    """Accumulate every variable while visiting each input file once."""
    reference_elements = _bulk_variable_elements(readers[0])
    cell_count = np.asarray(readers[0].read_variable("CellID")).size
    component_count = sum(
        int(reference_elements[name].attrib.get("vectorsize", 1))
        for name in selected
    )
    estimated_bytes = cell_count * component_count * np.dtype(np.float64).itemsize
    phase_start = time.perf_counter()
    progress_step = _progress_step(len(readers))
    print(
        "File-major float64 accumulators: approximately "
        f"{estimated_bytes / 1024**2:.1f} MiB "
        f"({estimated_bytes / 1024**3:.3f} GiB)",
        flush=True,
    )

    totals: dict[str, np.ndarray] = {}
    expected_shapes: dict[str, tuple[int, ...]] = {}

    for file_index, (reader, path, alignment) in enumerate(
        zip(readers, paths, alignments), start=1
    ):
        if _report_progress(file_index, len(readers), progress_step):
            print(
                f"  file-major read {file_index}/{len(readers)}: {path.name} "
                f"({time.perf_counter() - phase_start:.1f} s)",
                flush=True,
            )

        # Read arrays in physical byte-offset order for better file locality.
        file_elements = _bulk_variable_elements(reader)
        read_order = sorted(
            selected,
            key=lambda name: int(file_elements[name].text or 0),
        )
        for name in read_order:
            values = _read_aligned(reader, name, alignment)
            if file_index == 1:
                totals[name] = values.astype(np.float64, copy=True)
                expected_shapes[name] = values.shape
            else:
                if values.shape != expected_shapes[name]:
                    raise ValueError(
                        f"Shape of {name!r} differs in {path}: "
                        f"{values.shape} != {expected_shapes[name]}"
                    )
                np.add(totals[name], values, out=totals[name], casting="unsafe")

    for total in totals.values():
        total /= len(readers)
    return totals


def _initialize_average_worker(
    paths: Sequence[str], alignments: Sequence[np.ndarray | None]
) -> None:
    """Create process-owned readers and install read-only alignment arrays."""
    global _WORKER_READERS, _WORKER_ALIGNMENTS, _WORKER_PATHS
    _WORKER_PATHS = paths
    _WORKER_ALIGNMENTS = alignments
    _WORKER_READERS = [rhb.VlsvReader(path) for path in paths]


def _average_one_variable_worker(
    task: tuple[str, bool]
) -> tuple[str, np.ndarray, float]:
    """Average one variable in a worker, preserving serial file order."""
    name, report_file_progress = task
    if (
        _WORKER_READERS is None
        or _WORKER_ALIGNMENTS is None
        or _WORKER_PATHS is None
    ):
        raise RuntimeError("Multiprocessing worker was not initialized")

    variable_start = time.perf_counter()
    total = _read_aligned(
        _WORKER_READERS[0], name, _WORKER_ALIGNMENTS[0]
    ).astype(np.float64, copy=True)
    expected_shape = total.shape

    total_files = len(_WORKER_READERS)
    progress_step = _progress_step(total_files)
    if report_file_progress:
        print(
            f"  representative variable {name}: 1/{total_files} files",
            flush=True,
        )

    for file_index, (reader, path, alignment) in enumerate(
        zip(_WORKER_READERS[1:], _WORKER_PATHS[1:], _WORKER_ALIGNMENTS[1:]),
        start=2,
    ):
        values = _read_aligned(reader, name, alignment)
        if values.shape != expected_shape:
            raise ValueError(
                f"Shape of {name!r} differs in {path}: "
                f"{values.shape} != {expected_shape}"
            )
        np.add(total, values, out=total, casting="unsafe")

        if report_file_progress and _report_progress(
            file_index, total_files, progress_step
        ):
            print(
                f"  representative variable {name}: "
                f"{file_index}/{total_files} files "
                f"({time.perf_counter() - variable_start:.1f} s)",
                flush=True,
            )

    total /= len(_WORKER_READERS)
    return name, total, time.perf_counter() - variable_start


def _write_parallel_variable_major(
    readers: Sequence[rhb.VlsvReader],
    paths: Sequence[Path],
    output: Path,
    alignments: Sequence[np.ndarray | None],
    selected: Sequence[str],
    elements: dict[str, object],
    workers: int,
) -> None:
    """Average independent variables in processes and write in the parent."""
    worker_count = min(workers, len(selected))
    start_method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
    context = mp.get_context(start_method)
    path_strings = [str(path) for path in paths]
    representative = max(
        selected,
        key=lambda name: int(elements[name].attrib.get("vectorsize", 1)),
    )
    print(
        f"Multiprocessing: {worker_count} workers, {start_method} start method; "
        f"file-progress monitor {representative}",
        flush=True,
    )

    # Submit work before opening the output.  With fork, this prevents workers
    # from inheriting the writer's mutable file descriptor and XML state.
    with ProcessPoolExecutor(
        max_workers=worker_count,
        mp_context=context,
        initializer=_initialize_average_worker,
        initargs=(path_strings, alignments),
    ) as executor:
        tasks = [(name, name == representative) for name in selected]
        results = executor.map(_average_one_variable_worker, tasks, chunksize=1)
        with rhb.VlsvWriter(
            readers[0], str(output), copy_meshes=["SpatialGrid"]
        ) as writer:
            # executor.map preserves selected-variable order.  Each worker also
            # sums files in serial input order, preserving numerical results.
            for index, (name, mean, worker_elapsed) in enumerate(results, start=1):
                print(
                    f"[{index}/{len(selected)}] completed {name} in "
                    f"{worker_elapsed:.1f} s; writing",
                    flush=True,
                )
                _write_average(writer, elements, name, mean)


def _balanced_variable_groups(
    selected: Sequence[str], elements: dict[str, object], workers: int
) -> list[list[str]]:
    """Greedily balance groups by stored component count."""
    group_count = min(workers, len(selected))
    groups: list[list[str]] = [[] for _ in range(group_count)]
    weights = [0] * group_count
    order = sorted(
        selected,
        key=lambda name: (-int(elements[name].attrib.get("vectorsize", 1)), name),
    )
    for name in order:
        target = min(range(group_count), key=lambda index: (weights[index], index))
        groups[target].append(name)
        weights[target] += int(elements[name].attrib.get("vectorsize", 1))
    return groups


def _initialize_streaming_worker(
    paths: Sequence[str],
    reference_keys: np.ndarray,
    minimum_cellid: int,
    cell_count: int,
    chunk_cells: int,
) -> None:
    """Install streaming state and allocate reusable uint32 mapping arrays."""
    global _STREAM_PATHS
    global _STREAM_REFERENCE_KEYS
    global _STREAM_MINIMUM_CELLID
    global _STREAM_CELL_COUNT
    global _STREAM_CHUNK_CELLS
    global _STREAM_RAW_POSITIONS
    global _STREAM_POSITIONS

    _STREAM_PATHS = paths
    _STREAM_REFERENCE_KEYS = reference_keys
    _STREAM_MINIMUM_CELLID = minimum_cellid
    _STREAM_CELL_COUNT = cell_count
    _STREAM_CHUNK_CELLS = chunk_cells
    _STREAM_RAW_POSITIONS = np.arange(cell_count, dtype=np.uint32)
    _STREAM_POSITIONS = np.empty(cell_count, dtype=np.uint32)


def _streaming_alignment(
    reader: rhb.VlsvReader, path: str
) -> np.ndarray | None:
    """Build one dense uint32 alignment and discard it after this file."""
    if (
        _STREAM_REFERENCE_KEYS is None
        or _STREAM_MINIMUM_CELLID is None
        or _STREAM_CELL_COUNT is None
        or _STREAM_RAW_POSITIONS is None
        or _STREAM_POSITIONS is None
    ):
        raise RuntimeError("Streaming worker was not initialized")

    cellids = np.asarray(reader.read_variable("CellID"), dtype=np.int64)
    if cellids.ndim != 1 or cellids.size != _STREAM_CELL_COUNT:
        raise ValueError(f"CellID array size differs in {path}")
    cellids -= _STREAM_MINIMUM_CELLID
    if int(cellids.min()) < 0 or int(cellids.max()) >= _STREAM_CELL_COUNT:
        raise ValueError(f"CellID set differs in {path}")
    if np.array_equal(cellids, _STREAM_REFERENCE_KEYS):
        return None

    sentinel = np.iinfo(np.uint32).max
    _STREAM_POSITIONS.fill(sentinel)
    _STREAM_POSITIONS[cellids] = _STREAM_RAW_POSITIONS
    alignment = _STREAM_POSITIONS[_STREAM_REFERENCE_KEYS]
    if bool(np.any(alignment == sentinel)):
        raise ValueError(f"Duplicate or missing CellIDs in {path}")
    return alignment


def _streaming_update(
    accumulator: np.memmap,
    values: np.ndarray,
    alignment: np.ndarray | None,
    *,
    initialize: bool,
) -> None:
    """Copy or add one variable in chunks, avoiding a full aligned copy."""
    if _STREAM_CHUNK_CELLS is None:
        raise RuntimeError("Streaming worker was not initialized")
    for start in range(0, accumulator.shape[0], _STREAM_CHUNK_CELLS):
        stop = min(start + _STREAM_CHUNK_CELLS, accumulator.shape[0])
        source = (
            values[start:stop]
            if alignment is None
            else values[alignment[start:stop]]
        )
        target = accumulator[start:stop]
        if initialize:
            np.copyto(target, source, casting="unsafe")
        else:
            np.add(target, source, out=target, casting="unsafe")


def _streaming_group_worker(
    task: tuple[int, Sequence[tuple[str, str, tuple[int, ...]]], bool]
) -> tuple[int, list[str], float, int]:
    """Average a balanced variable group with per-file transient alignments."""
    if _STREAM_PATHS is None or _STREAM_CHUNK_CELLS is None:
        raise RuntimeError("Streaming worker was not initialized")
    group_index, specifications, report_progress = task
    started = time.perf_counter()
    accumulators = {
        name: np.memmap(path, dtype=np.float64, mode="r+", shape=shape)
        for name, path, shape in specifications
    }
    shapes = {name: shape for name, _, shape in specifications}
    names = [name for name, _, _ in specifications]
    progress_step = _progress_step(len(_STREAM_PATHS))
    reordered_files = 0

    for file_index, path in enumerate(_STREAM_PATHS, start=1):
        reader = rhb.VlsvReader(path)
        alignment = None if file_index == 1 else _streaming_alignment(reader, path)
        reordered_files += alignment is not None

        file_elements = _bulk_variable_elements(reader)
        read_order = sorted(
            names, key=lambda name: int(file_elements[name].text or 0)
        )
        for name in read_order:
            values = np.asarray(reader.read_variable(name))
            if values.shape != shapes[name]:
                raise ValueError(
                    f"Shape of {name!r} differs in {path}: "
                    f"{values.shape} != {shapes[name]}"
                )
            _streaming_update(
                accumulators[name],
                values,
                alignment,
                initialize=file_index == 1,
            )
            del values
        del alignment, reader

        if report_progress and _report_progress(
            file_index, len(_STREAM_PATHS), progress_step
        ):
            print(
                f"  streaming files {file_index}/{len(_STREAM_PATHS)}: "
                f"{Path(path).name} ({time.perf_counter() - started:.1f} s)",
                flush=True,
            )

    for accumulator in accumulators.values():
        for start in range(0, accumulator.shape[0], _STREAM_CHUNK_CELLS):
            stop = min(start + _STREAM_CHUNK_CELLS, accumulator.shape[0])
            target = accumulator[start:stop]
            np.divide(target, len(_STREAM_PATHS), out=target)
        accumulator.flush()
    del accumulators
    return group_index, names, time.perf_counter() - started, reordered_files


def _write_streaming_grouped(
    readers: Sequence[rhb.VlsvReader],
    paths: Sequence[Path],
    output: Path,
    reference_keys: np.ndarray,
    minimum_cellid: int,
    selected: Sequence[str],
    elements: dict[str, object],
    workers: int,
    chunk_cells: int,
    work_dir: os.PathLike[str] | str | None,
) -> None:
    """Run balanced streaming workers backed by temporary memory maps."""
    cell_count = int(reference_keys.size)
    groups = _balanced_variable_groups(selected, elements, workers)
    component_count = sum(
        int(elements[name].attrib.get("vectorsize", 1)) for name in selected
    )
    accumulator_bytes = cell_count * component_count * np.dtype(np.float64).itemsize
    mapping_bytes_per_worker = cell_count * (
        3 * np.dtype(np.uint32).itemsize + np.dtype(np.int64).itemsize
    )

    if work_dir is None:
        work_base = Path(
            os.environ.get("SLURM_TMPDIR") or tempfile.gettempdir()
        ).expanduser().resolve()
    else:
        work_base = Path(work_dir).expanduser().resolve()
    work_base.mkdir(parents=True, exist_ok=True)
    free_bytes = shutil.disk_usage(work_base).free
    if free_bytes < int(accumulator_bytes * 1.05):
        raise OSError(
            f"Insufficient free space in streaming work directory {work_base}: "
            f"need at least {_format_gib(int(accumulator_bytes * 1.05))}, "
            f"available {_format_gib(free_bytes)}"
        )

    weights = [
        sum(int(elements[name].attrib.get("vectorsize", 1)) for name in group)
        for group in groups
    ]
    print(
        f"Streaming strategy: {len(groups)} balanced groups with component "
        f"weights {weights}; accumulator workspace {_format_gib(accumulator_bytes)}",
        flush=True,
    )
    print(
        "Estimated peak CellID mapping storage: "
        f"{_format_gib(mapping_bytes_per_worker)} per worker, "
        f"{_format_gib(mapping_bytes_per_worker * len(groups))} aggregate",
        flush=True,
    )
    print(f"Streaming work directory: {work_base}", flush=True)
    print(f"Alignment/reorder chunk: {chunk_cells:,} cells", flush=True)

    start_method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
    context = mp.get_context(start_method)
    with tempfile.TemporaryDirectory(prefix="average_state_files_", dir=work_base) as temp:
        temp_path = Path(temp)
        specifications: dict[str, tuple[str, tuple[int, ...]]] = {}
        for index, name in enumerate(selected):
            vector_size = int(elements[name].attrib.get("vectorsize", 1))
            shape = (cell_count,) if vector_size == 1 else (cell_count, vector_size)
            accumulator_path = temp_path / f"{index:04d}.float64"
            mmap = np.memmap(
                accumulator_path, dtype=np.float64, mode="w+", shape=shape
            )
            mmap.flush()
            del mmap
            specifications[name] = (str(accumulator_path), shape)

        tasks = []
        for group_index, group in enumerate(groups):
            specs = [
                (name, specifications[name][0], specifications[name][1])
                for name in group
            ]
            tasks.append((group_index, specs, group_index == 0))

        path_strings = [str(path) for path in paths]
        with ProcessPoolExecutor(
            max_workers=len(groups),
            mp_context=context,
            initializer=_initialize_streaming_worker,
            initargs=(
                path_strings,
                reference_keys,
                minimum_cellid,
                cell_count,
                chunk_cells,
            ),
        ) as executor:
            results = executor.map(_streaming_group_worker, tasks, chunksize=1)
            for completed, (group_index, names, elapsed, reordered) in enumerate(
                results, start=1
            ):
                print(
                    f"[{completed}/{len(groups)}] streaming group {group_index + 1} "
                    f"complete in {elapsed:.1f} s; variables {', '.join(names)}; "
                    f"reordered files {reordered}/{len(paths)}",
                    flush=True,
                )

        print("Writing streaming accumulators to VLSV output", flush=True)
        with rhb.VlsvWriter(
            readers[0], str(output), copy_meshes=["SpatialGrid"]
        ) as writer:
            for index, name in enumerate(selected, start=1):
                accumulator_path, shape = specifications[name]
                mean = np.memmap(
                    accumulator_path, dtype=np.float64, mode="r", shape=shape
                )
                print(f"[{index}/{len(selected)}] writing {name}", flush=True)
                _write_average(writer, elements, name, mean)
                del mean


def average_state_files_files(
    input_files: Iterable[os.PathLike[str] | str],
    output_file: os.PathLike[str] | str,
    variables: Sequence[str] | None = None,
    *,
    overwrite: bool = False,
    strategy: str = "variable-major",
    workers: int = 1,
    chunk_cells: int = 1_000_000,
    work_dir: os.PathLike[str] | str | None = None,
) -> None:
    """Average selected bulk variables and write a new VLSV file.

    A float64 accumulator is used for every scalar or vector variable.  NaNs
    propagate according to ordinary NumPy arithmetic.  ``variable-major``
    minimizes memory use; ``file-major`` keeps every accumulator in memory to
    improve file locality when many or all variables are selected;
    ``streaming`` groups variables across workers and creates transient uint32
    CellID alignments backed by temporary memory-mapped accumulators.
    """
    if strategy not in {"variable-major", "file-major", "streaming"}:
        raise ValueError(
            "strategy must be 'variable-major', 'file-major', or 'streaming'"
        )
    if workers < 1:
        raise ValueError("workers must be at least 1")
    if chunk_cells < 1:
        raise ValueError("chunk_cells must be at least 1")
    if strategy == "file-major" and workers != 1:
        raise ValueError(
            "Multiprocessing is unavailable with file-major strategy; use "
            "--strategy variable-major, --strategy streaming, or --workers 1"
        )
    run_start = time.perf_counter()
    started_utc = datetime.now(timezone.utc)
    paths = [Path(path).expanduser().resolve() for path in input_files]
    if not paths:
        raise ValueError("At least one input VLSV file is required")
    if len(paths) != len(set(paths)):
        raise ValueError("Input VLSV file list contains duplicates")
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    output = Path(output_file).expanduser().resolve()
    log_path = Path(f"{output}.log")
    if output in paths:
        raise ValueError("The output file must not also be an input file")
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output already exists: {output}; use --force")
    if log_path.exists() and not overwrite:
        raise FileExistsError(f"Run log already exists: {log_path}; use --force")
    if overwrite:
        # A failed replacement run must not leave a successful-looking log
        # describing an older version of the output file.
        log_path.unlink(missing_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)

    print("RHybridPy VLSV averaging", flush=True)
    print(f"Input files: {len(paths)}", flush=True)
    print(f"Output file: {output}", flush=True)
    print(f"Strategy: {strategy}; requested workers: {workers}", flush=True)

    print("Phase 1/3: reading VLSV metadata", flush=True)
    metadata_start = time.perf_counter()
    metadata_step = _progress_step(len(paths))
    readers = []
    for file_index, path in enumerate(paths, start=1):
        readers.append(rhb.VlsvReader(str(path)))
        if _report_progress(file_index, len(paths), metadata_step):
            print(
                f"  metadata {file_index}/{len(paths)}: {path.name} "
                f"({time.perf_counter() - metadata_start:.1f} s)",
                flush=True,
            )
    print(
        f"Phase 1/3 complete in {time.perf_counter() - metadata_start:.1f} s",
        flush=True,
    )

    if strategy == "streaming":
        reference_keys, minimum_cellid, selected, elements = (
            _validate_streaming_inputs(readers, paths, variables)
        )
        alignments: list[np.ndarray | None] = []
    else:
        alignments, selected, elements = _validate_inputs(
            readers, paths, variables
        )
        reference_keys = np.empty(0, dtype=np.uint32)
        minimum_cellid = 0

    vector_sizes = [int(elements[name].attrib.get("vectorsize", 1)) for name in selected]
    scalar_count = sum(size == 1 for size in vector_sizes)
    vector_count = sum(size > 1 for size in vector_sizes)
    component_count = sum(vector_sizes)
    cell_count = (
        int(reference_keys.size)
        if strategy == "streaming"
        else np.asarray(readers[0].read_variable("CellID")).size
    )
    print(
        f"Grid cells: {cell_count:,}; selected variables: {len(selected)} "
        f"({scalar_count} scalar, {vector_count} vector, "
        f"{component_count} total components)",
        flush=True,
    )

    effective_workers = min(workers, len(selected))
    print(
        f"Phase 3/3: averaging and writing with {effective_workers} worker(s)",
        flush=True,
    )
    averaging_start = time.perf_counter()
    if strategy == "streaming":
        _write_streaming_grouped(
            readers,
            paths,
            output,
            reference_keys,
            minimum_cellid,
            selected,
            elements,
            effective_workers,
            chunk_cells,
            work_dir,
        )
    elif strategy == "variable-major" and effective_workers > 1:
        _write_parallel_variable_major(
            readers,
            paths,
            output,
            alignments,
            selected,
            elements,
            effective_workers,
        )
    else:
        if workers > 1:
            print("Only one variable selected; using one process", flush=True)

        # VlsvWriter opens with w+b, so an explicitly authorized overwrite is safe.
        with rhb.VlsvWriter(
            readers[0], str(output), copy_meshes=["SpatialGrid"]
        ) as writer:
            if strategy == "file-major":
                means = _average_file_major(readers, paths, alignments, selected)
                for index, name in enumerate(selected, start=1):
                    print(
                        f"[{index}/{len(selected)}] writing {name}", flush=True
                    )
                    _write_average(writer, elements, name, means[name])
            else:
                for index, name in enumerate(selected, start=1):
                    variable_start = time.perf_counter()
                    print(
                        f"[{index}/{len(selected)}] averaging {name}", flush=True
                    )

                    # The first file is already in reference/output order.
                    total = _read_aligned(
                        readers[0], name, alignments[0]
                    ).astype(np.float64, copy=True)
                    expected_shape = total.shape

                    # Each CellID permutation is reused for all variables.
                    for reader, path, alignment in zip(
                        readers[1:], paths[1:], alignments[1:]
                    ):
                        values = _read_aligned(reader, name, alignment)
                        if values.shape != expected_shape:
                            raise ValueError(
                                f"Shape of {name!r} differs in {path}: "
                                f"{values.shape} != {expected_shape}"
                            )
                        np.add(total, values, out=total, casting="unsafe")

                    total /= len(readers)
                    _write_average(writer, elements, name, total)
                    print(
                        f"[{index}/{len(selected)}] completed {name} in "
                        f"{time.perf_counter() - variable_start:.1f} s",
                        flush=True,
                    )

    print(
        f"Phase 3/3 complete in {time.perf_counter() - averaging_start:.1f} s",
        flush=True,
    )
    elapsed_seconds = time.perf_counter() - run_start
    run_log = _write_run_log(
        output=output,
        input_paths=paths,
        selected=selected,
        cell_count=cell_count,
        strategy=strategy,
        requested_workers=workers,
        effective_workers=effective_workers,
        chunk_cells=chunk_cells,
        started_utc=started_utc,
        elapsed_seconds=elapsed_seconds,
        overwrite=overwrite,
    )
    print(
        f"Wrote {len(selected)} averaged variables from {len(paths)} files "
        f"to {output} using {strategy} strategy and {effective_workers} worker(s); "
        f"total elapsed {elapsed_seconds:.1f} s",
        flush=True,
    )
    print(f"Wrote provenance log to {run_log}", flush=True)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_files",
        nargs="*",
        metavar="INPUT.vlsv",
        help=(
            "input VLSV files; the first file defines output CellID order; "
            "use --input-list for very large sets"
        ),
    )
    parser.add_argument(
        "--input-list",
        metavar="FILE",
        help="text file containing one input VLSV path per line",
    )
    parser.add_argument("-o", "--output", required=True, help="output VLSV file")
    parser.add_argument(
        "--variables",
        type=_parse_variables,
        help="comma-separated subset; default: all SpatialGrid bulk variables",
    )
    parser.add_argument(
        "--strategy",
        choices=("variable-major", "file-major", "streaming"),
        default="variable-major",
        help=(
            "averaging loop: variable-major minimizes memory (default); "
            "file-major keeps all float64 accumulators in RAM; streaming "
            "uses balanced workers, transient uint32 alignments, and temporary "
            "memory-mapped accumulators"
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        metavar="N",
        help=(
            "processes for variable-major or streaming averaging (default: 1); "
            "values above the number of selected variables are capped"
        ),
    )
    parser.add_argument(
        "--chunk-cells",
        type=int,
        default=1_000_000,
        metavar="N",
        help="cells per aligned update in streaming mode (default: 1000000)",
    )
    parser.add_argument(
        "--work-dir",
        help=(
            "temporary accumulator directory for streaming mode; default: "
            "$SLURM_TMPDIR when set, otherwise the system temporary directory"
        ),
    )
    parser.add_argument(
        "-f", "--force", action="store_true", help="overwrite an existing output"
    )
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()
    input_files = list(args.input_files)
    if args.input_list:
        try:
            input_files.extend(_read_input_list(args.input_list))
        except (OSError, ValueError) as error:
            parser.error(str(error))
    if not input_files:
        parser.error("provide INPUT.vlsv arguments or --input-list FILE")
    average_state_files_files(
        input_files,
        args.output,
        variables=args.variables,
        overwrite=args.force,
        strategy=args.strategy,
        workers=args.workers,
        chunk_cells=args.chunk_cells,
        work_dir=args.work_dir,
    )


if __name__ == "__main__":
    main()
