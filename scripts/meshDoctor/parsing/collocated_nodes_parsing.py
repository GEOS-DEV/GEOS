import logging
from typing import List, Dict

from checks import COLLOCATES_NODES
from checks.collocated_nodes import Options, Result

from . import cli_parsing

__TOLERANCE = "tolerance"


def get_help():
    return f"""Checks if nodes are collocated.
    options:
        {__TOLERANCE}: The absolute distance between two nodes for them to be considered collocated."""


def parse_options(options: Dict[str, str]) -> Options:
    """
    From the parsed cli options, return the converted options for collocated nodes check.
    :param options: Parsed cli options.
    :return: Options instance.
    """
    cli_parsing.validate_cli_options(COLLOCATES_NODES, {__TOLERANCE}, options)
    return Options(float(options[__TOLERANCE]))


def display_results(options: Options, result: Result):
    all_duplicated_nodes = []
    for bucket in result.nodes_buckets:
        for node in bucket:
            all_duplicated_nodes.append(node)
    all_duplicated_nodes = set(all_duplicated_nodes)  # Surely useless
    logging.warning(f"You have {len(all_duplicated_nodes)} collocated nodes (tolerance = {options.tolerance}).")

    logging.info(f"Here are all the buckets of collocated nodes.")
    tmp = []
    for bucket in result.nodes_buckets:
        tmp.append(f"({', '.join(map(str, bucket))})")
    logging.info(f"({', '.join(tmp)})")
