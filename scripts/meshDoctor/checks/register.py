from . import all_checks
from . import collocated_nodes, COLLOCATES_NODES


def register():
    """
    Register all the checks.
    :return:
    """
    all_checks[COLLOCATES_NODES] = collocated_nodes.check
