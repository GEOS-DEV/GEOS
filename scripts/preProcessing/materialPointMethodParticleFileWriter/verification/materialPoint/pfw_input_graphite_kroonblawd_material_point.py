#!/usr/bin/env python3
"""Generate a dry-run Graphite material-point case inspired by Kroonblawd Fig. 1.

Run this after building the optional C++ target by passing --run and the driver
executable path.  Without --run it writes the material XML and load-path CSV so
the generated input can be inspected.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

_THIS_DIR = Path(__file__).resolve().parents[2]
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from pfw_material_point import kroonblawd_graphite_case, run_material_point  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", default="geos_mpm_material_point_driver")
    parser.add_argument("--work-dir", default="graphite_kroonblawd_material_point")
    parser.add_argument("--theta-deg", type=float, default=75.0)
    parser.add_argument("--pressure-gpa", type=float, default=30.0)
    parser.add_argument("--dt", type=float, default=5.0e-7)
    parser.add_argument("--n-steps", type=int, default=1000)
    parser.add_argument("--strain-rate", type=float, default=-1.0e3)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()

    case = kroonblawd_graphite_case(
        theta_deg=args.theta_deg,
        pressure_gpa=args.pressure_gpa,
        dt=args.dt,
        n_steps=args.n_steps,
        strain_rate_us=args.strain_rate,
    )
    result = run_material_point(
        case,
        executable=args.driver,
        work_dir=args.work_dir,
        dry_run=not args.run,
    )
    print("material XML:", result.material_xml)
    print("load path CSV:", result.load_path_csv)
    print("output CSV:", result.output_csv)
    print("command:", " ".join(result.command))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
