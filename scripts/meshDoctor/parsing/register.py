import logging
from typing import Dict, Callable

from . import all_checks_helpers, CheckHelper

import checks.register
from checks import COLLOCATES_NODES


def __register_collocated_nodes() -> CheckHelper:
    from . import collocated_nodes_parsing
    return CheckHelper(parse_cli_options=collocated_nodes_parsing.parse_cli_options,
                       display_results=collocated_nodes_parsing.display_results,
                       get_help=collocated_nodes_parsing.get_help)


__HELPERS: Dict[str, Callable[[None], CheckHelper]] = {
    COLLOCATES_NODES: __register_collocated_nodes
}


def register():
    """
    Register all the parsing checks. Eventually initiate the registration of all the checks too.
    :return: None
    """
    loaded_checks = checks.register.register()
    for check_name in loaded_checks:
        all_checks_helpers[check_name] = __HELPERS[check_name]()
        logging.debug(f"Parsing for check \"{check_name}\" is loaded.")
