#!/usr/bin/env python3
"""Reject reference data members in objects copied into GEOS device kernels.

Copying an object containing a reference preserves the referenced address; it
does not move the referent to device memory (see GEOS PR #4194). Store scalar
parameters by value and LvArray views by value instead.

This source check covers classes/structs named KernelWrapper, names ending in
Update or Updates, and types named by KernelWrapper aliases. It also inspects
their known base classes and members stored by value. Constructor parameters,
method return types, local variables, and static members are outside this rule.

Tree-sitter parses source without compiling GEOS, including inactive #ifdef
branches. This is a syntax rule, not a proof of device accessibility: pointer
allocation, generated macro members, and template argument substitution are
not analyzed. Install dependencies with:
    python3 -m pip install -r scripts/code_checks_requirements.txt

With no arguments, check GEOS sources tracked by Git and untracked, nonignored
sources under src/coreComponents. Optional arguments select individual files.
"""

import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
import re
import subprocess
import sys

try:
    from tree_sitter import Language, Parser
    import tree_sitter_cpp
except ImportError:
    sys.exit("Missing C++ parser. Run: python3 -m pip install -r "
             "scripts/code_checks_requirements.txt")


SOURCE_SUFFIXES = (".hpp", ".cpp", ".hpp.template", ".cpp.template")
DECORATORS = re.compile(
    rb"\b(?:GEOS_HOST_DEVICE|GEOS_HOST|GEOS_DEVICE|GEOS_FORCE_INLINE)\b")


@dataclass
class Record:
    name: tuple
    path: Path
    bases: list = field(default_factory=list)
    members: list = field(default_factory=list)


def type_name(node):
    """Read a source type name, without template arguments."""
    if node is None:
        return None
    if node.type == "template_type":
        return type_name(node.child_by_field_name("name"))
    if node.type == "qualified_identifier":
        scope = node.child_by_field_name("scope")
        name = type_name(node.child_by_field_name("name"))
        if name is None:
            return None
        prefix = type_name(scope) if scope else ""
        if prefix is None:
            return None
        return prefix + "::" + name
    if node.type in ("type_identifier", "namespace_identifier", "primitive_type"):
        return node.text.decode()
    return None


def declarator_child(node):
    child = node.child_by_field_name("declarator")
    if child is None:
        child = next((item for item in node.named_children
                      if item.type.endswith("declarator") or
                      item.type in ("identifier", "field_identifier", "type_identifier")), None)
    return child


def binding_operator(node):
    """Find the innermost declarator operator (reference, pointer, function...).

    This distinguishes a reference member from a method returning a reference
    or a function pointer returning a reference. Parentheses do not change it.
    """
    operator = None
    while node is not None:
        if node.type.endswith("reference_declarator"):
            operator = "reference"
        elif node.type.endswith("pointer_declarator"):
            operator = "pointer"
        elif node.type.endswith("function_declarator"):
            operator = "function"
        elif node.type.endswith("array_declarator"):
            operator = "array"
        node = declarator_child(node)
    return operator


def lookup(name, scope, index):
    if not name:
        return None
    parts = tuple(name.lstrip(":").split("::"))
    prefixes = [()] if name.startswith("::") else [scope[:i] for i in range(len(scope), -1, -1)]
    for prefix in prefixes:
        key = prefix + parts
        if key in index:
            return key
    return None


def parse_sources(sources):
    parser = Parser(Language(tree_sitter_cpp.language()))
    records = defaultdict(list)
    aliases = defaultdict(list)

    for path, source in sources:
        # Preserve byte positions and line numbers when masking GEOS annotations.
        normalized = DECORATORS.sub(lambda match: b" " * len(match[0]), source)
        tree = parser.parse(normalized)

        def visit(node, scope=(), record=None):
            if node.type == "namespace_definition":
                name = node.child_by_field_name("name")
                if name:
                    scope += tuple(name.text.decode().split("::"))
            elif node.type in ("class_specifier", "struct_specifier"):
                name = type_name(node.child_by_field_name("name"))
                if name and node.child_by_field_name("body"):
                    scope += tuple(name.split("::"))
                    record = Record(scope, Path(path))
                    records[scope].append(record)
                    for child in node.named_children:
                        if child.type == "base_class_clause":
                            record.bases.extend(filter(None, (type_name(base) for base in child.named_children)))
            elif node.type == "alias_declaration":
                name = node.child_by_field_name("name")
                descriptor = node.child_by_field_name("type")
                if name and descriptor:
                    aliases[scope + (name.text.decode(),)].append((
                        type_name(descriptor.child_by_field_name("type")),
                        binding_operator(descriptor.child_by_field_name("declarator"))))
            elif node.type == "type_definition":
                for declarator in node.children_by_field_name("declarator"):
                    name = declarator
                    while declarator_child(name) is not None:
                        name = declarator_child(name)
                    if name.type == "type_identifier":
                        aliases[scope + (name.text.decode(),)].append((
                            type_name(node.child_by_field_name("type")), binding_operator(declarator)))
            elif node.type == "field_declaration" and record:
                is_static = any(child.type == "storage_class_specifier" and child.text == b"static"
                                for child in node.named_children)
                if not is_static:
                    for declarator in node.children_by_field_name("declarator"):
                        operator = binding_operator(declarator)
                        if operator != "function":
                            record.members.append((
                                type_name(node.child_by_field_name("type")), operator,
                                node.start_point.row + 1,
                                source[node.start_byte:node.end_byte].decode().strip()))
            for child in node.named_children:
                visit(child, scope, record)

        visit(tree.root_node)
    return records, aliases


def check_sources(sources):
    """Return (violations, number of inspected records) for (path, bytes) pairs."""
    records, aliases = parse_sources(sources)
    index = records.keys() | aliases.keys()

    def resolve(name, scope, seen=frozenset()):
        key = lookup(name, scope, index)
        if key is None or key in seen:
            return set(), False
        targets = {key} if key in records else set()
        is_reference = False
        for target, operator in aliases.get(key, []):
            if operator == "reference":
                is_reference = True
            elif operator not in ("pointer", "function"):
                nested, reference = resolve(target, key[:-1], seen | {key})
                targets |= nested
                is_reference |= reference
        return targets, is_reference

    pending = [key for key in index if key[-1] == "KernelWrapper" or
               key[-1].endswith(("Update", "Updates"))]
    checked = set()
    violations = set()
    while pending:
        key = pending.pop()
        if key in checked:
            continue
        checked.add(key)
        targets, _ = resolve("::" + "::".join(key), ())
        pending.extend(targets - checked)
        for record in records.get(key, []):
            for base in record.bases:
                targets, _ = resolve(base, key[:-1])
                pending.extend(targets - checked)
            for name, operator, line, declaration in record.members:
                if operator in ("pointer", "function"):
                    continue
                targets, reference = resolve(name, key)
                if operator == "reference" or reference:
                    violations.add((str(record.path), line, "::".join(key), declaration))
                else:
                    pending.extend(targets - checked)
    return sorted(violations), len(checked & records.keys())


def source_paths(repo_root):
    # Git lists a submodule as a gitlink; it does not descend into vendored code.
    result = subprocess.run(
        ["git", "-C", str(repo_root), "ls-files", "--cached", "--others",
         "--exclude-standard", "-z", "--", "src/coreComponents"],
        check=True, capture_output=True)
    return [repo_root / name.decode() for name in result.stdout.split(b"\0")
            if name and name.decode().endswith(SOURCE_SUFFIXES)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="*", type=Path, help="Source files to check instead of the GEOS tree")
    args = parser.parse_args()
    try:
        paths = args.files or source_paths(Path(__file__).resolve().parents[1])
        violations, count = check_sources((path, path.read_bytes()) for path in paths)
    except (OSError, subprocess.CalledProcessError) as error:
        print(f"Could not check device kernel members: {error}", file=sys.stderr)
        return 2
    for path, line, owner, declaration in violations:
        print(f"{path}:{line}: error: reference data member in device kernel type {owner}\n"
              f"  {' '.join(declaration.split())}\n"
              "  Store parameters or LvArray views by value; copying a reference preserves its address.")
    print(f"Device kernel member check: {len(violations)} violation(s), "
          f"{count} type(s) inspected, {len(paths)} file(s).")
    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main())
