import logging
from typing import Dict, Callable

from . import all_checks_helpers, CheckHelper, build_check_helper

import checks.register
from checks import COLLOCATES_NODES, GENERATE_FRACTURES, GENERATE_GLOBAL_IDS


def __register_collocated_nodes() -> CheckHelper:
    from . import collocated_nodes_parsing
    return build_check_helper(collocated_nodes_parsing)


def __register_generate_fractures() -> CheckHelper:
    from . import generate_fractures_parsing
    return build_check_helper(generate_fractures_parsing)


def __register_generate_global_ids() -> CheckHelper:
    from . import generate_global_ids_parsing
    return build_check_helper(generate_global_ids_parsing)


__HELPERS: Dict[str, Callable[[None], CheckHelper]] = {
    COLLOCATES_NODES: __register_collocated_nodes,
    GENERATE_FRACTURES: __register_generate_fractures,
    GENERATE_GLOBAL_IDS: __register_generate_global_ids,
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
