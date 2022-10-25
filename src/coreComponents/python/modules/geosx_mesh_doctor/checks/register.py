import logging
from typing import Tuple, Dict, Callable, Any

from . import all_checks, COLLOCATES_NODES, GENERATE_FRACTURES, GENERATE_GLOBAL_IDS


def __get_collocated_nodes_check():
    """
    Returns the "collocated nodes" check.
    :return: The "check" function.
    """
    from . import collocated_nodes
    return collocated_nodes.check


def __get_generate_fractures_check():
    """
    Returns the "collocated nodes" check.
    :return: The "check" function.
    """
    from . import generate_fractures
    return generate_fractures.check


def __get_generate_global_ids_check():
    """
    Returns the "collocated nodes" check.
    :return: The "check" function.
    """
    from . import generate_global_ids
    return generate_global_ids.check


__CHECKS: Dict[str, Callable[[None], Any]] = {
    COLLOCATES_NODES: __get_collocated_nodes_check,
    GENERATE_FRACTURES: __get_generate_fractures_check,
    GENERATE_GLOBAL_IDS: __get_generate_global_ids_check,
}


def register() -> Tuple[str]:
    """
    Register all the checks.
    A check that fails to load won't stop the process.
    :return: The keys of all the registered checks.
    """
    loaded_checks = []
    for check_name, check_provider in __CHECKS.items():
        try:
            all_checks[check_name] = check_provider()
            loaded_checks.append(check_name)
            logging.debug(f"Check \"{check_name}\" is loaded.")
        except Exception as e:
            logging.warning(f"Could not load module \"{check_name}\": {e}")
    return loaded_checks
