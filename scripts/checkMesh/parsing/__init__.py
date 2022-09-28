from dataclasses import dataclass
from typing import Callable, Dict, Any


@dataclass(frozen=True)
class CheckHelper:
    parse: Callable[[Dict[str, str]], Any]
    display_results: Callable[[Any, Any], None]
    get_help: Callable[[None], str]


all_checks_helpers = dict()
