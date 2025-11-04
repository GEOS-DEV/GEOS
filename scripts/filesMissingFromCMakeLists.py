#!/usr/bin/env python3
import os
import re
from pathlib import Path

# ----------------------------
# config
# ----------------------------

# Directories to skip when scanning the tree for cpp files
SKIP_DIR_NAMES = {
    '.git', '.github', '.gitlab', '.idea', '.vscode', '__pycache__',
    'build', 'build-debug', 'build-release', 'cmake-build-debug',
    'cmake-build-release', 'install', 'out', 'dist','blt'
}

# File extensions considered "source files"
SRC_EXTENSIONS = {'.cpp', '.hpp'}

# CMake file names to parse
CMAKE_FILE_NAMES = {'CMakeLists.txt', '.cmake'}  # second entry lets us also catch foo.cmake if you extend script


# ----------------------------
# helpers
# ----------------------------

def is_skipped_dir(dirname: str) -> bool:
    """Return True if we should not descend into this directory."""
    return dirname in SKIP_DIR_NAMES


def collect_all_cpp_files(repo_root: Path) -> set[Path]:
    """
    Walk the repo and return every *.cpp-like file path (relative to repo_root)
    except in skipped dirs.
    """
    all_cpp = set()

    for dirpath, dirnames, filenames in os.walk(repo_root):
        # prune skip dirs in-place so walk does not descend
        dirnames[:] = [d for d in dirnames if not is_skipped_dir(d)]

        for f in filenames:
            p = Path(dirpath) / f
            if p.suffix in SRC_EXTENSIONS:
                all_cpp.add(p.relative_to(repo_root))

    return all_cpp


def find_all_cmake_files(repo_root: Path) -> list[Path]:
    """
    Return list of all CMakeLists.txt and *.cmake files (except in skipped dirs).
    """
    cmake_files = []
    for dirpath, dirnames, filenames in os.walk(repo_root):
        dirnames[:] = [d for d in dirnames if not is_skipped_dir(d)]

        for f in filenames:
            if f == "CMakeLists.txt" or f.endswith(".cmake"):
                cmake_files.append(Path(dirpath) / f)

    return cmake_files


def tokenize_cmake_sources(cmake_text: str) -> list[str]:
    """
    Super light tokenizer:
    - strips comments (# ...)
    - splits on any whitespace, parens, quotes
    - returns tokens that look like file paths
    """
    # remove comments
    no_comments = []
    for line in cmake_text.splitlines():
        # CMake comments start with '#'
        if '#' in line:
            line = line.split('#', 1)[0]
        no_comments.append(line)
    text = "\n".join(no_comments)

    # We want to capture things that look like paths to source files.
    # We'll just extract stuff ending in .cpp/.cc/.cxx/.C etc using regex.
    pattern = r'([A-Za-z0-9_./\\+-]+(?:\.(?:cpp|cc|cxx|C)))'
    return re.findall(pattern, text)


def normalize_and_filter(tokens: list[str], cmake_dir: Path, repo_root: Path) -> set[Path]:
    """
    Convert token strings to normalized relative Paths if they look like real files.
    Handles:
      - relative paths like src/foo.cpp
      - absolute paths under repo root
    Ignores anything that doesn't exist.
    """
    out = set()
    for tok in tokens:
        raw = Path(tok)

        # If token is relative, interpret it relative to the CMake file location
        if not raw.is_absolute():
            cand = (cmake_dir / raw).resolve()
        else:
            cand = raw.resolve()

        try:
            rel = cand.relative_to(repo_root.resolve())
        except ValueError:
            # file is outside repo_root
            continue

        if cand.exists() and cand.suffix in SRC_EXTENSIONS:
            out.add(rel)

    return out


def main():
    repo_root = Path(os.getcwd()).resolve()

    # 1. gather cpp files from disk
    fs_cpp = collect_all_cpp_files(repo_root)

    # 2. gather cpp files from cmake
    cmake_files = find_all_cmake_files(repo_root)
    cmake_cpp: set[Path] = set()

    for cmake_path in cmake_files:
        try:
            txt = cmake_path.read_text()
        except Exception as e:
            print(f"Warning: could not read {cmake_path}: {e}")
            continue

        tokens = tokenize_cmake_sources(txt)
        cmake_cpp |= normalize_and_filter(tokens, cmake_path.parent, repo_root)

    # 3. diff
    unused_cpp = sorted(fs_cpp - cmake_cpp)
    missing_on_disk = sorted(cmake_cpp - fs_cpp)  # sanity check

    # 4. report
    print("=== Summary ===")
    print(f"Total source files on disk: {len(fs_cpp)}")
    print(f"Total source files referenced in CMake: {len(cmake_cpp)}")
    print()

    print("=== Present on disk but NOT referenced in any CMake file ===")
    if unused_cpp:
        for p in unused_cpp:
            print(p.as_posix())
    else:
        print("(none)")

    print()
    print("=== Referenced in CMake but file not found on disk (possible stale entry) ===")
    if missing_on_disk:
        for p in missing_on_disk:
            print(p.as_posix())
    else:
        print("(none)")


if __name__ == "__main__":
    main()
