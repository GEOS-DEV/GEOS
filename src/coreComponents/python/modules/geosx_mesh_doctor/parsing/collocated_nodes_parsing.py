import logging
import textwrap

from checks.collocated_nodes import Options, Result

from . import cli_parsing, COLLOCATES_NODES

__TOLERANCE = "tolerance"


def get_help():
    msg = f"""\
    Checks if nodes are collocated.
    
    {__TOLERANCE} [float]: The absolute distance between two nodes for them to be considered collocated.
    """
    return textwrap.dedent(msg)


def parse_cli_options(options_str: str) -> Options:
    """
    From the parsed cli options, return the converted options for collocated nodes check.
    :param options_str: Parsed cli options.
    :return: Options instance.
    """
    options = cli_parsing.parse_cli_option(options_str)
    cli_parsing.validate_cli_options(COLLOCATES_NODES, {__TOLERANCE}, options)
    return Options(float(options[__TOLERANCE]))


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
