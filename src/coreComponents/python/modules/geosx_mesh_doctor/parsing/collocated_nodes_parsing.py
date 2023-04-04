import logging

from checks.collocated_nodes import Options, Result

from . import COLLOCATES_NODES

__TOLERANCE = "tolerance"


def convert(parsed_options) -> Options:
    return Options(parsed_options[__TOLERANCE])


def fill_subparser(subparsers) -> None:
    p = subparsers.add_parser(COLLOCATES_NODES,
                              help="Checks if nodes are collocated.")
    p.add_argument('--' + __TOLERANCE,
                   type=float,
                   required=True,
                   help="[float]: The absolute distance between two nodes for them to be considered collocated.")


def display_results(options: Options, result: Result):
    all_duplicated_nodes = []
    for bucket in result.nodes_buckets:
        for node in bucket:
            all_duplicated_nodes.append(node)
    all_duplicated_nodes = set(all_duplicated_nodes)    # Surely useless
    if all_duplicated_nodes:
        logging.error(f"You have {len(all_duplicated_nodes)} collocated nodes (tolerance = {options.tolerance}).")

        logging.info("Here are all the buckets of collocated nodes.")
        tmp = []
        for bucket in result.nodes_buckets:
            tmp.append(f"({', '.join(map(str, bucket))})")
        logging.info(f"({', '.join(tmp)})")
    else:
        logging.error(f"You have no collocated node (tolerance = {options.tolerance}).")

    if result.wrong_support_elements:
        tmp = ", ".join(map(str, result.wrong_support_elements))
        logging.error(f"You have {len(result.wrong_support_elements)} elements with duplicated support nodes.\n" + tmp)
    else:
        logging.error("You have no element with duplicated support nodes.")
