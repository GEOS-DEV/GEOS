import logging
from typing import Dict, Callable, Any, Tuple

from parsing import (
    CheckHelper,
    COLLOCATES_NODES,
    ELEMENT_VOLUMES,
    FIX_ELEMENTS_ORDERINGS,
    GENERATE_FRACTURES,
    GENERATE_GLOBAL_IDS,
    NON_CONFORMAL,
    SELF_INTERSECTING_ELEMENTS,
    SUPPORTED_ELEMENTS,
)

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
    import parsing.element_volumes_parsing
    return __build_check_helper(parsing.element_volumes_parsing)


def __get_element_volumes_check():
    from checks import element_volumes
    return element_volumes.check


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


def __get_non_conformal_check_helper() -> CheckHelper:
    import parsing.non_conformal_parsing
    return __build_check_helper(parsing.non_conformal_parsing)


def __get_non_conformal_check():
    from checks import non_conformal
    return non_conformal.check


__HELPERS[NON_CONFORMAL] = __get_non_conformal_check_helper
__CHECKS[NON_CONFORMAL] = __get_non_conformal_check


def __get_self_intersecting_elements_check_helper() -> CheckHelper:
    import parsing.self_intersecting_elements_parsing
    return __build_check_helper(parsing.self_intersecting_elements_parsing)


def __get_self_intersecting_elements_check():
    from checks import self_intersecting_elements
    return self_intersecting_elements.check


__HELPERS[SELF_INTERSECTING_ELEMENTS] = __get_self_intersecting_elements_check_helper
__CHECKS[SELF_INTERSECTING_ELEMENTS] = __get_self_intersecting_elements_check


def __get_supported_elements_check_helper() -> CheckHelper:
    import parsing.supported_elements_parsing
    return __build_check_helper(parsing.supported_elements_parsing)


def __get_supported_elements_check():
    from checks import supported_elements
    return supported_elements.check


__HELPERS[SUPPORTED_ELEMENTS] = __get_supported_elements_check_helper
__CHECKS[SUPPORTED_ELEMENTS] = __get_supported_elements_check


def __get_fix_elements_orderings_helper() -> CheckHelper:
    import parsing.fix_elements_orderings_parsing
    return __build_check_helper(parsing.fix_elements_orderings_parsing)


def __get_fix_elements_orderings_check():
    from checks import fix_elements_orderings
    return fix_elements_orderings.check


__HELPERS[FIX_ELEMENTS_ORDERINGS] = __get_fix_elements_orderings_helper
__CHECKS[FIX_ELEMENTS_ORDERINGS] = __get_fix_elements_orderings_check


# TODO use importlib

def __load_checks() -> Dict[str, Callable[[str, Any], Any]]:
    """
    Loads all the checks.
    This function acts like a protection layer if a module fails to load.
    A check that fails to load won't stop the process.
    :return: The checks.
    """
    loaded_checks: Dict[str, Callable[[str, Any], Any]] = dict()
    for check_name, check_provider in __CHECKS.items():
        try:
            loaded_checks[check_name] = check_provider()
            logging.debug(f"Check \"{check_name}\" is loaded.")
        except Exception as e:
            logging.warning(f"Could not load module \"{check_name}\": {e}")
    return loaded_checks


def register() -> Tuple[Dict[str, Callable[[str, Any], Any]], Dict[str, CheckHelper]]:
    """
    Register all the parsing checks. Eventually initiate the registration of all the checks too.
    :return: The checks and the checks helpers.
    """
    loaded_checks: Dict[str, Callable[[str, Any], Any]] = __load_checks()
    loaded_checks_helpers: Dict[str, CheckHelper] = dict()
    for check_name in loaded_checks.keys():
        loaded_checks_helpers[check_name] = __HELPERS[check_name]()
        logging.debug(f"Parsing for check \"{check_name}\" is loaded.")
    return loaded_checks, loaded_checks_helpers
