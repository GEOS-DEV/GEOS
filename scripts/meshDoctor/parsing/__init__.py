from dataclasses import dataclass
from typing import Callable, Any


@dataclass(frozen=True)
class CheckHelper:
    parse_cli_options: Callable[[str], Any]
    display_results: Callable[[Any, Any], None]
    get_help: Callable[[None], str]


all_checks_helpers = dict()
