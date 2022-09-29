from . import all_checks_helpers, CheckHelper
from . import collocated_nodes_parsing

from checks import COLLOCATES_NODES
from checks import register as checks_register


def register():
    """
    Register all the parsing checks. Eventually initiate the registration of all the checks too.
    :return: None
    """
    all_checks_helpers[COLLOCATES_NODES] = CheckHelper(parse_cli_options=collocated_nodes_parsing.parse_cli_options,
                                                       display_results=collocated_nodes_parsing.display_results,
                                                       get_help=collocated_nodes_parsing.get_help)
    checks_register.register()  # TODO do not fail if a check failed to load.
