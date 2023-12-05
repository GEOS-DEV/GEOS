import argparse
from dataclasses import dataclass

from typing import (
    Iterator,
    Sequence,
)

import pytest

from checks.vtk_utils import (
    VtkOutput,
)

from checks.generate_fractures import (
    FracturePolicy,
    Options,
)
from parsing.generate_fractures_parsing import (
    convert,
    display_results,
    fill_subparser,
)


@dataclass(frozen=True)
class TestCase:
    __test__ = False
    cli_args: Sequence[str]
    options: Options
    exception: bool = False


def __generate_generate_fractures_parsing_test_data() -> Iterator[TestCase]:
    field: str = "attribute"
    main_mesh: str = "output.vtu"
    fracture_mesh: str = "fracture.vtu"

    cli_gen: str = f"generate_fractures --policy {{}} --name {field} --values 0,1 --output {main_mesh} --fracture-output {fracture_mesh}"

    cli_args: Sequence[str] = cli_gen.format("field").split()
    options: Options = Options(policy=FracturePolicy.FIELD, field=field, field_values=frozenset((0, 1)),
                               vtk_output=VtkOutput(output=main_mesh, is_data_mode_binary=True),
                               vtk_fracture_output=VtkOutput(output=fracture_mesh, is_data_mode_binary=True))
    yield TestCase(cli_args, options)

    cli_args: Sequence[str] = cli_gen.format("internal_surfaces").split()
    options: Options = Options(policy=FracturePolicy.INTERNAL_SURFACES, field=field, field_values=frozenset((0, 1)),
                               vtk_output=VtkOutput(output=main_mesh, is_data_mode_binary=True),
                               vtk_fracture_output=VtkOutput(output=fracture_mesh, is_data_mode_binary=True))
    yield TestCase(cli_args, options)
    cli_args: Sequence[str] = cli_gen.format("dummy").split()
    options: Options = Options(policy=FracturePolicy.INTERNAL_SURFACES, field=field, field_values=frozenset((0, 1)),
                               vtk_output=VtkOutput(output=main_mesh, is_data_mode_binary=True),
                               vtk_fracture_output=VtkOutput(output=fracture_mesh, is_data_mode_binary=True))
    yield TestCase(cli_args, options, True)


def __f(test_case: TestCase):
    parser = argparse.ArgumentParser(description='Testing.')
    subparsers = parser.add_subparsers()
    fill_subparser(subparsers)
    args = parser.parse_args(test_case.cli_args)
    options = convert(vars(args))
    assert options.policy == test_case.options.policy
    assert options.field == test_case.options.field
    assert options.field_values == test_case.options.field_values


def test_display_results():
    # Dummy test for code coverage only. Shame on me!
    display_results(None, None)


@pytest.mark.parametrize("test_case", __generate_generate_fractures_parsing_test_data())
def test(test_case: TestCase):
    if test_case.exception:
        with pytest.raises(SystemExit):
            __f(test_case)
    else:
        __f(test_case)
