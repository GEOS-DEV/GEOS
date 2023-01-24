import importlib
import logging
from typing import Dict, Callable, Any, Tuple

import parsing
from parsing import CheckHelper


__HELPERS: Dict[str, Callable[[None], CheckHelper]] = dict()
__CHECKS: Dict[str, Callable[[None], Any]] = dict()


def __load_module_check(module_name: str, check_fct="check"):
    module = importlib.import_module("checks." + module_name)
    return getattr(module, check_fct)


def __load_module_check_helper(module_name: str, parsing_fct_suffix="_parsing"):
    module = importlib.import_module("parsing." + module_name + parsing_fct_suffix)
    return CheckHelper(parse_cli_options=module.parse_cli_options,
                       display_results=module.display_results,
                       get_help=module.get_help)


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
    def closure_trick(cn: str):
        __HELPERS[check_name] = lambda: __load_module_check_helper(cn)
        __CHECKS[check_name] = lambda: __load_module_check(cn)
    # Register the modules to load here.
    for check_name in (parsing.COLLOCATES_NODES,
                       parsing.ELEMENT_VOLUMES,
                       parsing.FIX_ELEMENTS_ORDERINGS,
                       parsing.GENERATE_FRACTURES,
                       parsing.GENERATE_GLOBAL_IDS,
                       parsing.NON_CONFORMAL,
                       parsing.SELF_INTERSECTING_ELEMENTS,
                       parsing.SUPPORTED_ELEMENTS):
        closure_trick(check_name)
    loaded_checks: Dict[str, Callable[[str, Any], Any]] = __load_checks()
    loaded_checks_helpers: Dict[str, CheckHelper] = dict()
    for check_name in loaded_checks.keys():
        loaded_checks_helpers[check_name] = __HELPERS[check_name]()
        logging.debug(f"Parsing for check \"{check_name}\" is loaded.")
    return loaded_checks, loaded_checks_helpers
