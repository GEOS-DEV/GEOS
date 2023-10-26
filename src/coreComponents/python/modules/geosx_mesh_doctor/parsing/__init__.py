import argparse
from dataclasses import dataclass
from typing import Callable, Any


COLLOCATES_NODES = "collocated_nodes"
ELEMENT_VOLUMES = "element_volumes"
FIX_ELEMENTS_ORDERINGS = "fix_elements_orderings"
GENERATE_CUBE = "generate_cube"
GENERATE_FRACTURES = "generate_fractures_2"
GENERATE_GLOBAL_IDS = "generate_global_ids"
NON_CONFORMAL = "non_conformal"
SELF_INTERSECTING_ELEMENTS = "self_intersecting_elements"
SUPPORTED_ELEMENTS = "supported_elements"


@dataclass(frozen=True)
class CheckHelper:
    fill_subparser: Callable[[Any], argparse.ArgumentParser]
    convert: Callable[[Any], Any]
    display_results: Callable[[Any, Any], None]
