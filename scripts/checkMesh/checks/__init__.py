from . import collocated_nodes

COLLOCATES_NODES = "collocated_nodes"

# This is a poor attempt at registering products in the factory...
all_checks = {
    COLLOCATES_NODES: collocated_nodes.check
}
