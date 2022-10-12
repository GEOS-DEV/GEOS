from dataclasses import dataclass
from typing import Callable, Any


@dataclass(frozen=True)
class CheckHelper:
    parse_cli_options: Callable[[str], Any]
    display_results: Callable[[Any, Any], None]
    get_help: Callable[[None], str]


def build_check_helper(module) -> CheckHelper:
    """
    If a module has functions `parse_cli_options`, `display_results` and `get_help`,
    the `CheckHelper` built from those functions.
    :param module: Any module
    :return: The CheckHelper instance.
    """
    return CheckHelper(parse_cli_options=module.parse_cli_options,
                       display_results=module.display_results,
                       get_help=module.get_help)


all_checks_helpers = dict()
