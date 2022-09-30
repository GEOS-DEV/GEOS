import logging
from typing import Tuple, Dict, Callable, Any

from . import all_checks, COLLOCATES_NODES


def __get_collocated_nodes_check():
    """
    Returns the "collocated nodes" check.
    :return: The "check" function.
    """
    from . import collocated_nodes
    return collocated_nodes.check


__CHECKS: Dict[str, Callable[[None], Any]] = {
    COLLOCATES_NODES: __get_collocated_nodes_check
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
