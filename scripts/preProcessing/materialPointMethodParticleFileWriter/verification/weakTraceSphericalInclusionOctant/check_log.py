"""Check recorded grid-projection acceptance; this is not a stress-error test."""
import argparse
import json
import math
import re
from pathlib import Path


def parse_count_map(text):
    if text.strip() == "none":
        return {}
    counts = {}
    for item in text.split(","):
        if not item:
            continue
        key, value = item.split("=", 1)
        counts[key.strip()] = int(value)
    return counts


def sum_counts(counts):
    return sum(int(value) for value in counts.values())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log",type=Path)
    parser.add_argument("--maximum-jump",type=float,default=1e-10)
    parser.add_argument("--maximum-bc-residual",type=float,default=1e-11)
    parser.add_argument("--maximum-fallback-dofs",type=int)
    parser.add_argument("--maximum-requested-fallback-nodes",type=int)
    parser.add_argument("--maximum-fallback-node-fields",type=int)
    parser.add_argument("--minimum-time",type=float)
    args=parser.parse_args()
    text=args.log.read_text(errors="replace")
    rows=[]
    geometry_rows=[]
    for line in text.splitlines():
        if "C0 boundary-compatible geometry:" in line:
            def integer(key):
                match=re.search(re.escape(key)+r"=([0-9]+)",line)
                if match is None: raise ValueError("Missing geometry log field: "+key)
                return int(match.group(1))
            request_match=re.search(r"fallback requests=\{([^}]*)\}",line)
            node_field_match=re.search(r"fallback node-fields=\{([^}]*)\}",line)
            requests=parse_count_map(request_match.group(1)) if request_match else {}
            node_fields=parse_count_map(node_field_match.group(1)) if node_field_match else {}
            geometry_rows.append(dict(
                candidate_rows=integer("candidate rows"),
                requested_fallback_nodes=integer("requested fallback nodes"),
                fallback_requests=requests,
                fallback_node_fields=node_fields,
                fallback_request_total=sum_counts(requests),
                fallback_node_field_total=sum_counts(node_fields),
            ))
        if "C0 boundary-compatible QR component=" not in line:
            continue
        def number(key):
            match=re.search(re.escape(key)+r"=([^,\s]+)",line)
            if match is None: raise ValueError("Missing log field: "+key)
            value=float(match.group(1))
            if not math.isfinite(value): raise ValueError("Nonfinite log field: "+key)
            return value
        rows.append(dict(component=int(number("component")),time=number("time"),
                         maximum_jump=number("max average jump"),boundary_residual=number("boundary residual"),
                         relative_residual=number("relative residual"),fallback_dofs=int(number("fallback dofs"))))
    if not rows:
        raise SystemExit("No compatible-QR records. Enable the MPMSubroutines log category/solver log level and confirm the mode is enabled.")
    if "C0 boundary-compatible projection failed" in text or "post-mortar boundary/synchronization update changed" in text:
        raise SystemExit("The log contains a failed coupled projection/update.")
    maximum_jump=max(r["maximum_jump"] for r in rows)
    maximum_bc=max(r["boundary_residual"] for r in rows)
    fallback=max(r["fallback_dofs"] for r in rows)
    last=max(r["time"] for r in rows)
    if maximum_jump>args.maximum_jump or maximum_bc>args.maximum_bc_residual:
        raise SystemExit("Grid constraint residual exceeds the requested diagnostic threshold.")
    if args.maximum_fallback_dofs is not None and fallback>args.maximum_fallback_dofs:
        raise SystemExit("Fallback footprint exceeds the requested threshold.")
    if args.maximum_requested_fallback_nodes is not None and geometry_rows:
        requested=max(r["requested_fallback_nodes"] for r in geometry_rows)
        if requested>args.maximum_requested_fallback_nodes:
            raise SystemExit("Requested fallback nodes exceed the requested threshold.")
    if args.maximum_fallback_node_fields is not None and geometry_rows:
        node_fields=max(r["fallback_node_field_total"] for r in geometry_rows)
        if node_fields>args.maximum_fallback_node_fields:
            raise SystemExit("Fallback node-fields exceed the requested threshold.")
    if args.minimum_time is not None and last<args.minimum_time-1e-8*max(1.,abs(args.minimum_time)):
        raise SystemExit("Recorded projections do not reach the requested final time.")
    result=dict(records=len(rows),last_projection_time=last,max_average_jump=maximum_jump,
        max_boundary_residual=maximum_bc,max_fallback_dofs=fallback,
        note="Grid projection only. Check particle motion, force balance and stress error separately.")
    if geometry_rows:
        result["geometry_records"]=len(geometry_rows)
        result["max_requested_fallback_nodes"]=max(r["requested_fallback_nodes"] for r in geometry_rows)
        result["max_fallback_request_total"]=max(r["fallback_request_total"] for r in geometry_rows)
        result["max_fallback_node_field_total"]=max(r["fallback_node_field_total"] for r in geometry_rows)
        result["last_geometry_record"]=geometry_rows[-1]
    print(json.dumps(result,indent=2))


if __name__=="__main__": main()
