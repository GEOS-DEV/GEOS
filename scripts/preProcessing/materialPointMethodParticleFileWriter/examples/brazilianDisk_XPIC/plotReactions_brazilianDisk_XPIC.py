#!/usr/bin/env python3
"""Plot Brazilian disk y-face reactions from reactionHistory.csv."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Plot +y/-y reactions vs time")
    parser.add_argument("--run-dir", default=".")
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--case-name", default="brazilianDisk_XPIC")
    return parser.parse_args()


def as_float(value):
    try:
        return float(value)
    except Exception:
        return None


def read_reactions(path: Path):
    with path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        has_header = any(token in sample.splitlines()[0].lower() for token in ("time", "rx", "ry"))
        if has_header:
            reader = csv.DictReader(handle)
            names = list(reader.fieldnames or [])
            lookup = {name.strip().lower(): name for name in names}
            time_key = lookup.get("time") or names[0]
            rym_key = lookup.get("ry-") or lookup.get("rym") or lookup.get("y-")
            ryp_key = lookup.get("ry+") or lookup.get("ryp") or lookup.get("y+")
            if rym_key is None or ryp_key is None:
                # Known PFW order: time,F00,F11,F22,Lx,Ly,Lz,Rx-,Rx+,Ry-,Ry+,Rz-,Rz+,L00,L11,L22
                rows = []
                for row in reader:
                    vals = [as_float(row.get(name, "")) for name in names]
                    if len(vals) >= 11:
                        rows.append((vals[0], vals[9], vals[10]))
                return rows
            rows = []
            for row in reader:
                rows.append((as_float(row[time_key]), as_float(row[rym_key]), as_float(row[ryp_key])))
            return rows
        reader = csv.reader(handle)
        rows = []
        for raw in reader:
            vals = [as_float(v) for v in raw]
            if len(vals) >= 11 and vals[0] is not None:
                rows.append((vals[0], vals[9], vals[10]))
        return rows


def main():
    args = parse_args()
    run_dir = Path(args.run_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    source = run_dir / "reactionHistory.csv"
    if not source.is_file():
        raise SystemExit(f"missing reaction history: {source}")
    rows = [(t, rym, ryp) for t, rym, ryp in read_reactions(source) if t is not None and rym is not None and ryp is not None]
    if not rows:
        raise SystemExit(f"no usable y-face reaction rows in {source}")
    t = [r[0] for r in rows]
    rym = [r[1] for r in rows]
    ryp = [r[2] for r in rows]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(t, rym, label="Ry-")
    ax.plot(t, ryp, label="Ry+")
    ax.set_xlabel("time")
    ax.set_ylabel("reaction")
    ax.set_title(f"{args.case_name}: y-face reactions")
    ax.grid(True, alpha=0.35)
    ax.legend(loc="best")
    fig.tight_layout()

    png = output_dir / f"{args.case_name}_reactions_y.png"
    csv_out = output_dir / f"{args.case_name}_reactions_y.csv"
    fig.savefig(png, dpi=180)
    with csv_out.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", "Ry-", "Ry+"])
        writer.writerows(rows)
    print(f"wrote {png}")
    print(f"wrote {csv_out}")


if __name__ == "__main__":
    main()
