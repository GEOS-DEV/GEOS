#!/usr/bin/env python3
"""Plot compression-direction reactions and assemble four ParticleColor frames."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default=".")
    parser.add_argument("--frames-dir", default=None)
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--case-name", default="rigidBodyJamming2D")
    parser.add_argument("--output-prefix", default="rigidBodyJamming2D")
    parser.add_argument("--compression-direction", choices=("x", "y", "z"), default="y")
    return parser.parse_args()


def as_float(value):
    try:
        return float(value)
    except Exception:
        return None


def compact(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())


def read_reaction_history(path: Path, direction: str):
    if not path.is_file():
        raise FileNotFoundError(f"missing reaction history: {path}")
    axis = direction.lower()
    minus_candidates = [f"r{axis}-", f"r{axis}m", f"{axis}-", f"reaction{axis}-", f"reaction{axis}minus"]
    plus_candidates = [f"r{axis}+", f"r{axis}p", f"{axis}+", f"reaction{axis}+", f"reaction{axis}plus"]
    fallback = {"x": (7, 8), "y": (9, 10), "z": (11, 12)}[axis]
    rows = []
    with path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        first = sample.splitlines()[0] if sample.splitlines() else ""
        has_header = any(token in first.lower() for token in ("time", "rx", "ry", "rz"))
        if has_header:
            reader = csv.DictReader(handle)
            names = list(reader.fieldnames or [])
            lookup = {compact(name): name for name in names}
            time_key = lookup.get("time") or (names[0] if names else None)
            minus_key = next((lookup.get(compact(c)) for c in minus_candidates if lookup.get(compact(c))), None)
            plus_key = next((lookup.get(compact(c)) for c in plus_candidates if lookup.get(compact(c))), None)
            for row in reader:
                if minus_key is None or plus_key is None:
                    values = [as_float(row.get(name, "")) for name in names]
                    if len(values) > max(fallback):
                        rows.append((values[0], values[fallback[0]], values[fallback[1]]))
                else:
                    rows.append((as_float(row.get(time_key, "")), as_float(row.get(minus_key, "")), as_float(row.get(plus_key, ""))))
        else:
            reader = csv.reader(handle)
            for raw in reader:
                values = [as_float(v) for v in raw]
                if len(values) > max(fallback):
                    rows.append((values[0], values[fallback[0]], values[fallback[1]]))
    return [r for r in rows if all(v is not None for v in r)]


def write_reaction_plot(run_dir: Path, output_dir: Path, prefix: str, direction: str):
    rows = read_reaction_history(run_dir / "reactionHistory.csv", direction)
    if not rows:
        raise RuntimeError("reactionHistory.csv did not contain usable rows")
    t = [r[0] for r in rows]
    rminus = [r[1] for r in rows]
    rplus = [r[2] for r in rows]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(t, rminus, label=f"R{direction}-")
    ax.plot(t, rplus, label=f"R{direction}+")
    ax.set_xlabel("time")
    ax.set_ylabel(f"{direction}-direction reaction")
    ax.set_title(f"{prefix}: compression-direction reactions")
    ax.grid(True, alpha=0.35)
    ax.legend(loc="best")
    fig.tight_layout()
    png = output_dir / f"{prefix}_reactions_{direction}.png"
    csv_out = output_dir / f"{prefix}_reactions_{direction}.csv"
    fig.savefig(png, dpi=180)
    with csv_out.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", f"R{direction}-", f"R{direction}+"])
        writer.writerows(rows)
    print(f"wrote {png}")
    print(f"wrote {csv_out}")


def make_color_sheet(frames_dir: Path, output_dir: Path, prefix: str):
    frames = sorted(frames_dir.glob(f"{prefix}_*_particleColor*.png"))
    labels = ["initial", "quarter", "threequarter", "final"]
    selected = []
    for label in labels:
        match = [p for p in frames if f"_{label}_" in p.name]
        if match:
            selected.append(match[-1])
    if len(selected) < 4:
        print("not enough ParticleColor frames for a four-time contact sheet; found " + str([p.name for p in frames]))
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.0))
    for ax, label, path in zip(axes.ravel(), labels, selected):
        ax.imshow(mpimg.imread(path))
        ax.set_title(label)
        ax.axis("off")
    fig.suptitle(f"{prefix}: ParticleColor at four times")
    fig.tight_layout()
    png = output_dir / f"{prefix}_color_four_times.png"
    fig.savefig(png, dpi=180)
    print(f"wrote {png}")


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).expanduser().resolve()
    frames_dir = Path(args.frames_dir).expanduser().resolve() if args.frames_dir else run_dir / "visit_frames"
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_reaction_plot(run_dir, output_dir, args.output_prefix, args.compression_direction)
    make_color_sheet(frames_dir, output_dir, args.output_prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
