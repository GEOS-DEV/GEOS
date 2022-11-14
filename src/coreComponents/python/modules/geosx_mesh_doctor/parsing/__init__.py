from dataclasses import dataclass
from typing import Dict, Callable, Any


COLLOCATES_NODES = "collocated_nodes"
ELEMENT_VOLUMES = "element_volumes"
GENERATE_FRACTURES = "generate_fractures"
GENERATE_GLOBAL_IDS = "generate_global_ids"


@dataclass(frozen=True)
class CheckHelper:
    parse_cli_options: Callable[[str], Any]
    display_results: Callable[[Any, Any], None]
    get_help: Callable[[None], str]


# Singleton-like pattern for check helpers registration
all_checks_helpers: Dict[str, CheckHelper] = dict()
