import logging
from typing import Dict, Callable, Any, Tuple

from parsing import (
    all_checks_helpers,
    CheckHelper,
    COLLOCATES_NODES,
    ELEMENT_VOLUMES,
    GENERATE_FRACTURES,
    GENERATE_GLOBAL_IDS,
)
from checks import all_checks

__HELPERS: Dict[str, Callable[[None], CheckHelper]] = dict()
__CHECKS: Dict[str, Callable[[None], Any]] = dict()


def __build_check_helper(module) -> CheckHelper:
    """
    If a module has functions `parse_cli_options`, `display_results` and `get_help`,
    the `CheckHelper` built from those functions.
    :param module: Any module
    :return: The CheckHelper instance.
    """
    return CheckHelper(parse_cli_options=module.parse_cli_options,
                       display_results=module.display_results,
                       get_help=module.get_help)


def __get_collocated_nodes_check_helper() -> CheckHelper:
    import parsing.collocated_nodes_parsing
    return __build_check_helper(parsing.collocated_nodes_parsing)


def __get_collocated_nodes_check():
    from checks import collocated_nodes
    return collocated_nodes.check


__HELPERS[COLLOCATES_NODES] = __get_collocated_nodes_check_helper
__CHECKS[COLLOCATES_NODES] = __get_collocated_nodes_check


def __get_element_volumes_check_helper() -> CheckHelper:
    import parsing.elements_volumes_parsing
    return __build_check_helper(parsing.elements_volumes_parsing)


def __get_element_volumes_check():
    from checks import elements_volumes
    return elements_volumes.check


__HELPERS[ELEMENT_VOLUMES] = __get_element_volumes_check_helper
__CHECKS[ELEMENT_VOLUMES] = __get_element_volumes_check


def __get_generate_fractures_check_helper() -> CheckHelper:
    import parsing.generate_fractures_parsing
    return __build_check_helper(parsing.generate_fractures_parsing)


def __get_generate_fractures_check():
    from checks import generate_fractures
    return generate_fractures.check


__HELPERS[GENERATE_FRACTURES] = __get_generate_fractures_check_helper
__CHECKS[GENERATE_FRACTURES] = __get_generate_fractures_check


def __get_generate_global_ids_check_helper() -> CheckHelper:
    import parsing.generate_global_ids_parsing
    return __build_check_helper(parsing.generate_global_ids_parsing)


def __get_generate_global_ids_check():
    from checks import generate_global_ids
    return generate_global_ids.check


__HELPERS[GENERATE_GLOBAL_IDS] = __get_generate_global_ids_check_helper
__CHECKS[GENERATE_GLOBAL_IDS] = __get_generate_global_ids_check


def __load_checks() -> Tuple[str]:
    """
    Loads all the checks.
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


def register():
    """
    Register all the parsing checks. Eventually initiate the registration of all the checks too.
    Fills the global variables `parsing.all_checks_helpers` and `checks.all_checks` with the loaded check helpers.
    :return: None
    """
    loaded_checks = __load_checks()
    for check_name in loaded_checks:
        all_checks_helpers[check_name] = __HELPERS[check_name]()
        logging.debug(f"Parsing for check \"{check_name}\" is loaded.")
