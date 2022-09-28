from dataclasses import dataclass
from typing import Callable, Dict, Any

from checks import COLLOCATES_NODES

from . import collocated_nodes_parsing

@dataclass
class CheckHelper:
    parse: Callable[[Dict[str, str]], Any]
    present: Callable[[Any, Any], None]
    get_help: Callable[[None], str]


all_checks_helpers = {
    COLLOCATES_NODES: CheckHelper(parse=collocated_nodes_parsing.parse_options,
                                  present=collocated_nodes_parsing.present_results,
                                  get_help=collocated_nodes_parsing.get_help)
}
