#!/usr/bin/env python3
"""Compare ATS result classifications from two integrated-test runs.

The ATS ``[Info]`` section contains run-specific timestamps and host names and
is intentionally ignored. Every field in the remaining sections is compared
as a semicolon-separated collection, so ordering differences from concurrent
ATS scheduling do not hide a test moving between result classifications.
"""

import argparse
import configparser
import sys
import unittest


def read_results(path):
    parser = configparser.ConfigParser(interpolation=None)
    parser.optionxform = str
    with open(path, encoding="utf-8") as stream:
        parser.read_file(stream)
    return parser


def normalize(value):
    return tuple(sorted(item.strip() for item in value.split(";") if item.strip()))


def differences(first, second):
    diffs = []
    sections = sorted(set(first.sections()) | set(second.sections()))
    for section in sections:
        if section.lower() == "info":
            continue
        first_options = set(first.options(section)) if first.has_section(section) else set()
        second_options = set(second.options(section)) if second.has_section(section) else set()
        for option in sorted(first_options | second_options):
            first_value = (first.get(section, option, raw=True)
                           if option in first_options else None)
            second_value = (second.get(section, option, raw=True)
                            if option in second_options else None)
            first_normalized = normalize(first_value) if first_value is not None else None
            second_normalized = normalize(second_value) if second_value is not None else None
            if first_normalized != second_normalized:
                diffs.append((section, option, first_normalized, second_normalized))
    return diffs


def compare(first_path, second_path):
    first = read_results(first_path)
    second = read_results(second_path)
    diffs = differences(first, second)
    if not diffs:
        print("ATS result classifications match (run metadata ignored).")
        return 0

    print(f"ATS result classifications differ in {len(diffs)} field(s):")
    for section, option, first_value, second_value in diffs:
        first_only = sorted(set(first_value or ()) - set(second_value or ()))
        second_only = sorted(set(second_value or ()) - set(first_value or ()))
        print(f"  [{section}] {option}: "
              f"first-only={len(first_only)}, second-only={len(second_only)}")
        if first_only:
            print(f"    first-only: {';'.join(first_only[:20])}")
        if second_only:
            print(f"    second-only: {';'.join(second_only[:20])}")
    return 1


class CompareIntegratedTestResultsTests(unittest.TestCase):
    def testNormalizeIgnoresOrderAndEmptyEntries(self):
        self.assertEqual(normalize("b;; a; b "), ("a", "b", "b"))

    def testInfoIsIgnored(self):
        first = configparser.ConfigParser(interpolation=None)
        second = configparser.ConfigParser(interpolation=None)
        first.read_dict({"Info": {"time": "one"}, "Results": {"passed": "a;b"}})
        second.read_dict({"Info": {"time": "two"}, "Results": {"passed": "b;a"}})
        self.assertEqual(differences(first, second), [])

    def testResultChangeIsDetected(self):
        first = configparser.ConfigParser(interpolation=None)
        second = configparser.ConfigParser(interpolation=None)
        first.read_dict({"Results": {"passed": "a", "failed": "b"}})
        second.read_dict({"Results": {"passed": "a;b", "failed": ""}})
        self.assertEqual(len(differences(first, second)), 2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("first", nargs="?")
    parser.add_argument("second", nargs="?")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        result = unittest.TextTestRunner(verbosity=2).run(
            unittest.defaultTestLoader.loadTestsFromTestCase(
                CompareIntegratedTestResultsTests))
        return 0 if result.wasSuccessful() else 1
    if not args.first or not args.second:
        parser.error("first and second ATS test_results.ini files are required")
    return compare(args.first, args.second)


if __name__ == "__main__":
    sys.exit(main())
