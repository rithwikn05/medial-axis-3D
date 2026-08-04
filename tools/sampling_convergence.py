#!/usr/bin/env python3
"""Run repeatable surface-sampling convergence experiments."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import time
from pathlib import Path


FIELDS = [
    "sampling_mode",
    "requested_samples",
    "actual_samples",
    "tetrahedra",
    "interior_voronoi_vertices",
    "contained_voronoi_edges",
    "rejected_exterior_edges",
    "raw_medial_polygons",
    "raw_medial_triangles",
    "unique_inward_poles",
    "validated_poles",
    "minimum_h_over_lfs",
    "maximum_h_over_lfs",
    "retained_candidate_polygons",
    "source_candidate_polygons",
    "supported_triangles",
    "supported_components",
    "final_components",
    "final_triangles",
    "retained_medial_area",
    "retained_boundary_length",
    "retained_seam_length",
    "retained_artificial_boundary_loops",
    "retained_artificial_boundary_length",
    "boundary_edges",
    "seam_edges",
    "junction_vertices",
    "boundary_loops",
    "artificial_boundary_edges",
    "elapsed_seconds",
    "return_code",
]


def first_match(text: str, pattern: str, groups: tuple[str, ...]) -> dict[str, str]:
    match = re.search(pattern, text, flags=re.MULTILINE)
    if match is None:
        return {group: "" for group in groups}
    return dict(zip(groups, match.groups()))


def parse_metrics(output: str) -> dict[str, str]:
    metrics: dict[str, str] = {}
    metrics.update(
        first_match(
            output,
            r"Loaded (\d+) samples, generated (\d+) tetrahedra, "
            r"\d+ boundary triangles, (\d+) interior Voronoi vertices, "
            r"(\d+) contained Voronoi edges \(rejected (\d+) "
            r"exterior-crossing edges\), and (\d+) medial-sheet polygons "
            r"\((\d+) triangles\)\.",
            (
                "actual_samples",
                "tetrahedra",
                "interior_voronoi_vertices",
                "contained_voronoi_edges",
                "rejected_exterior_edges",
                "raw_medial_polygons",
                "raw_medial_triangles",
            ),
        )
    )
    metrics.update(
        first_match(
            output,
            r"Selected (\d+) unique inward poles; (\d+) passed "
            r"medial-ball validation\.",
            ("unique_inward_poles", "validated_poles"),
        )
    )
    metrics.update(
        first_match(
            output,
            r"h/LFS density \[([-+0-9.eE]+), ([-+0-9.eE]+)\]",
            ("minimum_h_over_lfs", "maximum_h_over_lfs"),
        )
    )
    metrics.update(
        first_match(
            output,
            r"Retained (\d+) of (\d+) candidate Voronoi polygons as "
            r"medial sheets \((\d+) triangles in (\d+) components?\)",
            (
                "retained_candidate_polygons",
                "source_candidate_polygons",
                "supported_triangles",
                "supported_components",
            ),
        )
    )
    metrics.update(
        first_match(
            output,
            r"Component filtering retained (\d+) components? and (\d+) "
            r"triangles; removed \d+ triangles\. The retained complex has "
            r"(\d+) boundary edges, (\d+) non-manifold seam edges, and "
            r"(\d+) junction vertices, (\d+) detected boundary loops, and "
            r"(\d+) unresolved artificial boundary edges\.",
            (
                "final_components",
                "final_triangles",
                "boundary_edges",
                "seam_edges",
                "junction_vertices",
                "boundary_loops",
                "artificial_boundary_edges",
            ),
        )
    )
    component_pattern = re.compile(
        r"  Component \d+: \d+ triangles, area ([-+0-9.eE]+) "
        r"\([-+0-9.eE]+%\).*?boundary length ([-+0-9.eE]+), "
        r"seam length ([-+0-9.eE]+), boundary loops \d+ "
        r"\((\d+) artificial, ([-+0-9.eE]+)% of boundary length\) "
        r"-> (retained|removed)\.",
        flags=re.MULTILINE,
    )
    retained_area = 0.0
    retained_boundary_length = 0.0
    retained_seam_length = 0.0
    retained_artificial_loops = 0
    retained_artificial_boundary_length = 0.0
    for match in component_pattern.finditer(output):
        if match.group(6) != "retained":
            continue
        area = float(match.group(1))
        boundary_length = float(match.group(2))
        seam_length = float(match.group(3))
        artificial_loops = int(match.group(4))
        artificial_fraction = float(match.group(5)) / 100.0
        retained_area += area
        retained_boundary_length += boundary_length
        retained_seam_length += seam_length
        retained_artificial_loops += artificial_loops
        retained_artificial_boundary_length += (
            artificial_fraction * boundary_length
        )
    metrics["retained_medial_area"] = f"{retained_area:.17g}"
    metrics["retained_boundary_length"] = (
        f"{retained_boundary_length:.17g}"
    )
    metrics["retained_seam_length"] = f"{retained_seam_length:.17g}"
    metrics["retained_artificial_boundary_loops"] = str(
        retained_artificial_loops
    )
    metrics["retained_artificial_boundary_length"] = (
        f"{retained_artificial_boundary_length:.17g}"
    )
    return metrics


def default_viewer(repo_root: Path) -> Path:
    candidates = [
        repo_root / "build" / "Release" / "polyscope_viewer.exe",
        repo_root / "build" / "polyscope_viewer.exe",
        repo_root / "build" / "polyscope_viewer",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run the medial-axis viewer headlessly at several sample counts "
            "and write convergence metrics to CSV."
        )
    )
    parser.add_argument("input", type=Path, help="input TetGen .node file")
    parser.add_argument(
        "--samples",
        nargs="+",
        type=int,
        required=True,
        help="sample counts to test, such as: --samples 4500 6000 9000",
    )
    parser.add_argument(
        "--viewer",
        type=Path,
        help="polyscope_viewer executable (auto-detected by default)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="output CSV path (defaults under build/)",
    )
    parser.add_argument(
        "--min-contact-angle",
        type=float,
        help="optional contact-angle override passed to the viewer",
    )
    parser.add_argument(
        "--adaptive-sampling",
        action="store_true",
        help="pass --adaptive-sampling to the viewer",
    )
    parser.add_argument(
        "--profile-stages",
        action="store_true",
        help="pass --profile-stages to the viewer and preserve timings in logs",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        help="optional timeout in seconds for each run",
    )
    arguments = parser.parse_args()

    if any(count < 4 for count in arguments.samples):
        parser.error("every sample count must be at least 4")

    repo_root = Path(__file__).resolve().parent.parent
    input_path = arguments.input.resolve()
    viewer_path = (
        arguments.viewer.resolve()
        if arguments.viewer is not None
        else default_viewer(repo_root)
    )
    output_path = (
        arguments.output.resolve()
        if arguments.output is not None
        else repo_root
        / "build"
        / f"{input_path.stem}_sampling_convergence.csv"
    )
    log_directory = (
        output_path.parent / f"{output_path.stem}_logs"
    )

    if not input_path.is_file():
        parser.error(f"input does not exist: {input_path}")
    if not viewer_path.is_file():
        parser.error(f"viewer does not exist: {viewer_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_directory.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int | float]] = []
    for count in arguments.samples:
        command = [
            str(viewer_path),
            str(input_path),
            "--samples",
            str(count),
            "--no-cross-resolution",
            "--no-gui",
        ]
        if arguments.min_contact_angle is not None:
            command.extend([
                "--min-contact-angle",
                str(arguments.min_contact_angle),
            ])
        if arguments.adaptive_sampling:
            command.append("--adaptive-sampling")
        if arguments.profile_stages:
            command.append("--profile-stages")

        print(f"Running {count} requested samples...", flush=True)
        started = time.perf_counter()
        try:
            completed = subprocess.run(
                command,
                cwd=repo_root,
                capture_output=True,
                text=True,
                timeout=arguments.timeout,
                check=False,
            )
            elapsed = time.perf_counter() - started
            combined_output = completed.stdout
            if completed.stderr:
                combined_output += "\n[stderr]\n" + completed.stderr
            return_code = completed.returncode
        except subprocess.TimeoutExpired as error:
            elapsed = time.perf_counter() - started
            stdout = error.stdout or ""
            stderr = error.stderr or ""
            if isinstance(stdout, bytes):
                stdout = stdout.decode(errors="replace")
            if isinstance(stderr, bytes):
                stderr = stderr.decode(errors="replace")
            combined_output = stdout + "\n[stderr]\n" + stderr
            return_code = 124

        log_path = log_directory / f"samples_{count}.log"
        log_path.write_text(combined_output, encoding="utf-8")

        row: dict[str, str | int | float] = {
            field: "" for field in FIELDS
        }
        row.update(parse_metrics(combined_output))
        row["sampling_mode"] = (
            "adaptive" if arguments.adaptive_sampling else "uniform"
        )
        row["requested_samples"] = count
        row["elapsed_seconds"] = f"{elapsed:.6f}"
        row["return_code"] = return_code
        rows.append(row)

        with output_path.open("w", newline="", encoding="utf-8") as output:
            writer = csv.DictWriter(output, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)

        status = "ok" if return_code == 0 else f"exit {return_code}"
        print(f"  {status} in {elapsed:.2f}s; log: {log_path}", flush=True)

    print(f"Wrote convergence table: {output_path}")
    return 0 if all(row["return_code"] == 0 for row in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
