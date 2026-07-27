#!/usr/bin/env python3
"""Convert an ASCII OFF surface mesh to TetGen .node and .face files."""

from __future__ import annotations

import argparse
from pathlib import Path


def data_lines(path: Path):
    with path.open("r", encoding="utf-8-sig") as source:
        for line_number, raw_line in enumerate(source, 1):
            line = raw_line.split("#", 1)[0].strip()
            if line:
                yield line_number, line


def read_off(path: Path):
    lines = iter(data_lines(path))
    try:
        magic_line, magic = next(lines)
    except StopIteration as exc:
        raise ValueError("the input file is empty") from exc

    if magic == "OFF":
        try:
            count_line, count_text = next(lines)
        except StopIteration as exc:
            raise ValueError("the OFF counts line is missing") from exc
    elif magic.startswith("OFF "):
        count_line = magic_line
        count_text = magic[4:].strip()
    else:
        raise ValueError(f"line {magic_line}: expected the OFF header")

    counts = count_text.split()
    if len(counts) < 2:
        raise ValueError(
            f"line {count_line}: expected vertex and face counts"
        )

    try:
        vertex_count = int(counts[0])
        face_count = int(counts[1])
    except ValueError as exc:
        raise ValueError(f"line {count_line}: invalid OFF counts") from exc

    vertices: list[tuple[float, float, float]] = []
    for vertex_index in range(vertex_count):
        try:
            line_number, line = next(lines)
        except StopIteration as exc:
            raise ValueError(
                f"expected {vertex_count} vertices, found {vertex_index}"
            ) from exc
        fields = line.split()
        if len(fields) < 3:
            raise ValueError(
                f"line {line_number}: a vertex requires three coordinates"
            )
        try:
            vertices.append(
                (float(fields[0]), float(fields[1]), float(fields[2]))
            )
        except ValueError as exc:
            raise ValueError(
                f"line {line_number}: invalid vertex coordinates"
            ) from exc

    triangles: list[tuple[int, int, int]] = []
    for face_index in range(face_count):
        try:
            line_number, line = next(lines)
        except StopIteration as exc:
            raise ValueError(
                f"expected {face_count} faces, found {face_index}"
            ) from exc
        fields = line.split()
        try:
            corner_count = int(fields[0])
            corners = [int(value) for value in fields[1 : corner_count + 1]]
        except (IndexError, ValueError) as exc:
            raise ValueError(f"line {line_number}: invalid face") from exc

        if corner_count < 3 or len(corners) != corner_count:
            raise ValueError(
                f"line {line_number}: a face requires at least three vertices"
            )
        if any(index < 0 or index >= vertex_count for index in corners):
            raise ValueError(
                f"line {line_number}: face index is outside the vertex range"
            )

        for corner in range(1, corner_count - 1):
            triangle = (corners[0], corners[corner], corners[corner + 1])
            if len(set(triangle)) != 3:
                raise ValueError(
                    f"line {line_number}: face contains a degenerate triangle"
                )
            triangles.append(triangle)

    return vertices, triangles


def write_tetgen(
    output_base: Path,
    vertices: list[tuple[float, float, float]],
    triangles: list[tuple[int, int, int]],
) -> tuple[Path, Path]:
    node_path = output_base.with_suffix(".node")
    face_path = output_base.with_suffix(".face")
    node_path.parent.mkdir(parents=True, exist_ok=True)

    with node_path.open("w", encoding="utf-8", newline="\n") as output:
        output.write(f"{len(vertices)} 3 0 0\n")
        for index, (x, y, z) in enumerate(vertices, 1):
            output.write(f"{index} {x:.17g} {y:.17g} {z:.17g}\n")

    with face_path.open("w", encoding="utf-8", newline="\n") as output:
        output.write(f"{len(triangles)} 0\n")
        for index, (a, b, c) in enumerate(triangles, 1):
            output.write(f"{index} {a + 1} {b + 1} {c + 1}\n")

    return node_path, face_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert an ASCII OFF mesh to TetGen .node and .face files."
    )
    parser.add_argument("input", type=Path, help="input ASCII .off file")
    parser.add_argument(
        "output_base",
        nargs="?",
        type=Path,
        help="output path without an extension (defaults to the input path)",
    )
    arguments = parser.parse_args()

    input_path = arguments.input.resolve()
    output_base = (
        arguments.output_base.resolve()
        if arguments.output_base
        else input_path.with_suffix("")
    )

    try:
        vertices, triangles = read_off(input_path)
        node_path, face_path = write_tetgen(
            output_base, vertices, triangles
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))

    print(
        f"Wrote {len(vertices)} vertices and {len(triangles)} triangles:\n"
        f"  {node_path}\n"
        f"  {face_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
