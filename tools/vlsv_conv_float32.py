#!/usr/bin/env python3
"""Convert float64 physical-variable arrays into a compact float32 VLSV file.

Requirements
------------
Python 3.8+ and NumPy. RHybridPy is not required.
Only little-endian VLSV files with a zero endian marker are supported.

Usage
-----
Inspect the planned conversions and exact planned output size:
    python vlsv_conv_float32.py input.vlsv output.vlsv --dry-run

Convert using the default 64 MiB input chunks:
    python vlsv_conv_float32.py input.vlsv output.vlsv

Use smaller chunks and collect numerical rounding-error statistics:
    python vlsv_conv_float32.py input.vlsv output.vlsv --chunk-mib 16 --error-stats

Store only selected variables (CellID is always retained):
    python vlsv_conv_float32.py input.vlsv output.vlsv --variables cellB,cellE

Preview selection and output size without writing:
    python vlsv_conv_float32.py input.vlsv output.vlsv --variables cellB --dry-run

Explicitly allow replacement of an existing output:
    python vlsv_conv_float32.py input.vlsv output.vlsv --force

Behavior
--------
Converts XML VARIABLE entries with datatype="float" and datasize="8", except
CellID. By default all variables are retained. --variables accepts a comma-
separated list of exact, case-sensitive names (quote names containing spaces).
Unknown names and empty list entries are errors. Matching names on all meshes
are retained; selection does not distinguish meshes. CellID is always retained.
Unselected VARIABLE entries and their payloads are omitted, including integer
and existing float32 variables. All non-VARIABLE coordinate/mesh entries and
parameters, plus retained integer and float32 arrays, are copied byte-for-byte. XML attributes are preserved except converted
datasize values; all payload offsets are rebuilt.

The input is never modified and must not change during the run. The output
directory must exist and have space for the complete new file. With --force,
the old output also occupies space until the new file is ready. Conversion
writes a temporary file beside the destination, then publishes it on success.
The default no-overwrite publication requires filesystem hard-link support.

Float32 rounding is expected. Finite values overflowing to infinity cause
failure; nonzero values rounded to zero are counted. Existing NaNs and
infinities are accepted. --error-stats reports maximum absolute and relative
errors over finite converted values; relative error excludes original zeros.
Absolute errors across variables with different units are not comparable.

A 64 MiB input chunk uses 96 MiB for the input/output numeric buffers, plus
boolean checks, Python/XML overhead and filesystem caching. --error-stats uses
additional chunk-sized temporary arrays. Memory does not scale with cell count
apart from the small XML metadata. MiB means 1,048,576 bytes.

Progress throughput counts input payload bytes (not combined disk reads and
writes); ETA is approximate. Dry runs validate structure but do not scan values
for overflow. Exit status is 0 on success, 1 on handled conversion errors,
2 for command-line parsing errors, and 130 for a keyboard interrupt.
Temporary files are cleaned on handled failure; a forced kill may leave one.
"""

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import struct
import tempfile
import time
import xml.etree.ElementTree as ET

import numpy as np


@dataclass
class Block:
    """Original payload location and its associated mutable XML entry."""

    element: ET.Element
    offset: int
    size: int
    convert: bool

    @property
    def label(self):
        return f"{self.element.tag} {self.element.get('name', '')}"


def read_layout(stream):
    """Validate all payload bounds before creating output."""
    size = os.fstat(stream.fileno()).st_size
    # Bytes 0..7: endian marker; bytes 8..15: unsigned XML footer offset.
    header = stream.read(16)
    if len(header) != 16:
        raise ValueError('Truncated VLSV header')
    if header[:8] != bytes(8):
        raise ValueError('Unsupported endian marker; only little-endian VLSV is supported')
    footer_offset = struct.unpack('<Q', header[8:])[0]
    if not 16 <= footer_offset < size:
        raise ValueError('Invalid XML footer offset')
    stream.seek(footer_offset)
    root = ET.fromstring(stream.read())
    if root.tag != 'VLSV':
        raise ValueError('Expected a VLSV XML root')
    blocks = []
    for element in root:
        try:
            offset = int(element.text or '')
            count = int(element.attrib['arraysize'])
            components = int(element.attrib['vectorsize'])
            itemsize = int(element.attrib['datasize'])
        except (ValueError, KeyError) as error:
            raise ValueError(f'Invalid array metadata for {element.tag}') from error
        if count < 0 or components < 1 or itemsize < 1:
            raise ValueError(f'Invalid dimensions for {element.tag}')
        # arraysize counts rows; vectorsize counts scalar elements per row.
        length = count * components * itemsize
        if offset < 16 or offset + length > footer_offset:
            raise ValueError(f'Data block outside payload: {element.tag}')
        # Select by VLSV tag/type, never by the scientific variable name.
        convert = (element.tag == 'VARIABLE'
                   and element.get('name', '').lower() != 'cellid'
                   and element.get('datatype') == 'float' and itemsize == 8)
        blocks.append(Block(element, offset, length, convert))
    # Physical order keeps reads sequential even if XML entries are unordered.
    blocks.sort(key=lambda block: block.offset)
    end = 16
    for block in blocks:
        if block.size:
            if block.offset < end:
                raise ValueError(f'Overlapping data blocks: {block.label}')
            end = block.offset + block.size
    return root, blocks, size


def select_variables(root, blocks, variables):
    """Remove unselected VARIABLE entries, preserving CellID and mesh metadata.

    Selection applies to storage, independently of float64 conversion. Validate
    names before mutating the XML or creating any temporary output file.
    """
    if variables is None:
        return blocks, []
    names = [name.strip() for name in variables.split(',')]
    if any(not name for name in names):
        raise ValueError('--variables requires nonempty comma-separated names')
    requested = set(names)
    available = {block.element.get('name', '') for block in blocks
                 if block.element.tag == 'VARIABLE'}
    missing = requested - available
    if missing:
        raise ValueError('Unknown variable name(s): ' + ', '.join(sorted(missing))
                         + '. Available: ' + ', '.join(sorted(available)))
    kept, omitted = [], []
    for block in blocks:
        element = block.element
        name = element.get('name', '')
        if (element.tag == 'VARIABLE' and name.lower() != 'cellid'
                and name not in requested):
            root.remove(element)
            omitted.append(block)
        else:
            kept.append(block)
    return kept, omitted


def plan_footer(root, blocks):
    """Assign packed output offsets while retaining original Block locations."""
    # Element counts and vector dimensions stay unchanged; only byte sizes shrink.
    offset = 16
    for block in blocks:
        block.element.text = str(offset)
        if block.convert:
            block.element.set('datasize', '4')
        offset += block.size // 2 if block.convert else block.size
    return offset, ET.tostring(root, encoding='utf-8')


def convert_file(input_file, output_file, *, chunk_mib=64, dry_run=False,
                 force=False, error_stats=False, variables=None):
    """Validate, stream and publish a converted file; dry-run writes nothing."""
    source = Path(input_file).expanduser().resolve()
    destination = Path(output_file).expanduser().absolute()
    if not source.is_file():
        raise ValueError(f'Input is not a file: {source}')
    # Check path aliases, including symlinks and existing hard links.
    if destination.resolve() == source or (
        destination.exists() and os.path.samefile(source, destination)
    ):
        raise ValueError('Input and output must be different files')
    if os.path.lexists(destination) and not force:
        raise FileExistsError(f'Output already exists: {destination}; use --force')
    if not destination.parent.is_dir():
        raise ValueError(f'Output directory does not exist: {destination.parent}')
    if chunk_mib < 1:
        raise ValueError('chunk-mib must be at least 1')
    chunk_bytes = chunk_mib * 1024**2
    started = time.monotonic()
    with source.open('rb') as infile:
        initial_stat = os.fstat(infile.fileno())
        root, blocks, original_size = read_layout(infile)
        # Validate the complete input first, then omit unwanted payloads and XML.
        blocks, omitted = select_variables(root, blocks, variables)
        # Payload lengths determine every new offset without reading full arrays.
        footer_offset, footer = plan_footer(root, blocks)
        expected_size = footer_offset + len(footer)
        converted = sum(block.convert for block in blocks)
        print(f'Arrays: {converted} converted, {len(blocks)-converted} preserved, '
              f'{len(omitted)} omitted', flush=True)
        for block in omitted:
            print(f'  omit    {block.label}: {block.size:,} input bytes', flush=True)
        print(f'Size: {original_size:,} -> {expected_size:,} bytes '
              f'({100*(1-expected_size/original_size):.2f}% saved)', flush=True)
        for block in blocks:
            print(f"  {'convert' if block.convert else 'copy   '} {block.label}: "
                  f'{block.size:,} input bytes', flush=True)
        if dry_run:
            print('Dry run: no output written; numerical checks require conversion.')
            return
        total = sum(block.size for block in blocks)
        processed = underflow = 0
        last_report = started
        max_absolute = max_relative = 0.0
        temporary = None
        try:
            # Same-directory staging permits atomic publication on this filesystem.
            with tempfile.NamedTemporaryFile(mode='w+b', dir=destination.parent,
                                             prefix=f'.{destination.name}.',
                                             suffix='.tmp', delete=False) as outfile:
                temporary = Path(outfile.name)
                outfile.write(bytes(16))
                for index, block in enumerate(blocks, 1):
                    print(f'[{index}/{len(blocks)}] {block.label}', flush=True)
                    infile.seek(block.offset)
                    remaining = block.size
                    while remaining:
                        amount = min(remaining, chunk_bytes)
                        if block.convert:
                            # Flattened scalar chunks preserve cell/vector storage order.
                            values = np.fromfile(infile, dtype='<f8', count=amount//8)
                            if values.nbytes != amount:
                                raise ValueError(f'Truncated data: {block.label}')
                            with np.errstate(over='ignore', under='ignore', invalid='ignore'):
                                result = values.astype('<f4')
                            # Reject newly created infinities, but allow existing ones.
                            finite = np.isfinite(values)
                            if np.any(finite & np.isinf(result)):
                                raise ValueError(f'Float32 overflow in {block.label}')
                            underflow += int(np.count_nonzero((values != 0) & (result == 0)))
                            # Optional statistics cost extra memory and array passes.
                            if error_stats and np.any(finite):
                                original = values[finite]
                                difference = np.abs(result[finite].astype(np.float64)-original)
                                max_absolute = max(max_absolute, float(difference.max()))
                                nonzero = original != 0
                                if np.any(nonzero):
                                    max_relative = max(max_relative, float(
                                        (difference[nonzero]/np.abs(original[nonzero])).max()))
                                del original, difference, nonzero
                            result.tofile(outfile)
                            del values, result, finite
                        else:
                            # Raw copying preserves precision, datatype and exact bytes.
                            data = infile.read(amount)
                            if len(data) != amount:
                                raise ValueError(f'Truncated data: {block.label}')
                            outfile.write(data)
                            del data
                        remaining -= amount
                        processed += amount
                        now = time.monotonic()
                        if now-last_report >= 5 or processed == total:
                            elapsed = now-started
                            rate = processed/max(elapsed, 1e-9)
                            print(f'  {processed:,}/{total:,} input bytes; '
                                  f'{rate/1024**2:.1f} MiB/s; {elapsed:.1f} s elapsed; '
                                  f'ETA {(total-processed)/max(rate, 1e-9):.1f} s', flush=True)
                            last_report = now
                if outfile.tell() != footer_offset:
                    raise ValueError('Internal output-offset mismatch')
                # Write the footer once, then fill its offset into the header.
                outfile.write(footer)
                outfile.seek(8)
                outfile.write(struct.pack('<Q', footer_offset))
                outfile.flush()
                os.fsync(outfile.fileno())
            # Detect common input modifications before making the result visible.
            final_stat = os.fstat(infile.fileno())
            if (initial_stat.st_size, initial_stat.st_mtime_ns) != (
                final_stat.st_size, final_stat.st_mtime_ns
            ):
                raise ValueError('Input changed during conversion')
            if force:
                os.replace(temporary, destination)
            else:
                # Atomic publication without overwriting a concurrently created file.
                os.link(temporary, destination)
                temporary.unlink()
            temporary = None
        finally:
            # Keep an existing destination intact if conversion did not complete.
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    print(f'Nonzero values rounded to zero: {underflow:,}')
    if error_stats:
        print(f'Maximum absolute error: {max_absolute:.8g}; '
              f'maximum relative error: {max_relative:.8g} '
              '(finite converted values; absolute errors may have different units)')
    print(f'Completed: {destination}; {time.monotonic()-started:.1f} s; '
          f'{original_size:,} -> {expected_size:,} bytes; '
          f'{converted} converted, {len(blocks)-converted} preserved, '
          f'{len(omitted)} omitted arrays', flush=True)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('input_file', type=Path, help='source VLSV file (read only)')
    parser.add_argument('output_file', type=Path, help='new compact VLSV file')
    parser.add_argument('--variables', metavar='NAME[,NAME...]',
                        help='store only these exact variable names; always retain CellID '
                             'and non-variable data (default: retain all variables)')
    parser.add_argument('--chunk-mib', type=int, default=64,
                        help='input buffer size in MiB (default: 64)')
    parser.add_argument('--dry-run', action='store_true',
                        help='validate layout and show planned conversions without writing')
    parser.add_argument('--force', action='store_true', help='replace existing output')
    parser.add_argument('--error-stats', action='store_true',
                        help='measure rounding errors; uses extra memory and CPU')
    args = parser.parse_args()
    try:
        convert_file(**vars(args))
    except (OSError, ValueError, ET.ParseError) as error:
        parser.exit(1, f'ERROR: {error}\n')
    except KeyboardInterrupt:
        parser.exit(130, 'Interrupted; temporary output removed.\n')


if __name__ == '__main__':
    main()
