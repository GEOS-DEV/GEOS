#!/usr/bin/env python3
"""Extract a macOS TPL image archive into the given directory."""

import inspect
import os
import sys
import tarfile

def relocate(name):
    name = name[2:] if name.startswith("./") else name
    if name in ("", "."):
        return None
    if name == "spack-generated.cmake":
        return name
    prefix = "opt/GEOS/"
    if name.startswith(prefix):
        return name[len(prefix):]
    return None

root = sys.argv[2]
with tarfile.open(sys.argv[1]) as tar:
    members = tar.getmembers()
    kwargs = {}
    if "filter" in inspect.signature(tar.extractall).parameters:
        kwargs["filter"] = "fully_trusted"
    files = []
    for member in members:
        if member.islnk():
            continue
        relocated = relocate(member.name)
        if relocated is None:
            continue
        member.name = relocated
        files.append(member)
    try:
        tar.extractall(root, members=files, **kwargs)
    except OSError as exc:
        sys.stderr.write("EXTRACT_ERROR %s: %s\n" % (type(exc).__name__, exc))
        sys.exit(1)
    for member in members:
        if not member.islnk():
            continue
        relocated = relocate(member.name)
        link = relocate(member.linkname)
        if relocated is None or link is None:
            continue
        destination = os.path.join(root, relocated)
        source = os.path.join(root, link)
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        if os.path.lexists(destination):
            os.unlink(destination)
        os.link(source, destination)
