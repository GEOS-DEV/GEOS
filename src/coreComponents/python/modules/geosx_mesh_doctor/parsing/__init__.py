from dataclasses import dataclass
from typing import Callable, Any


COLLOCATES_NODES = "collocated_nodes"
ELEMENT_VOLUMES = "element_volumes"
GENERATE_FRACTURES = "generate_fractures"
GENERATE_GLOBAL_IDS = "generate_global_ids"
SELF_INTERSECTING_ELEMENTS = "self_intersecting_elements"


@dataclass(frozen=True)
class CheckHelper:
    parse_cli_options: Callable[[str], Any]
    display_results: Callable[[Any, Any], None]
    get_help: Callable[[None], str]
