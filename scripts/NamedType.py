import logging

import gdb

logger = logging.getLogger("NamedType")


class NamedTypePrinter(object):
    def __init__(self, val: gdb.Value):
        self.val = val["value_"]

    def to_string(self) -> str:
        return self.val


def build_array_printer():
    pp = gdb.printing.RegexpCollectionPrettyPrinter("NamedType")
    pp.add_printer('fluent::NamedType', '^fluent::NamedType<.*>$', NamedTypePrinter)
    return pp


try:
    import gdb.printing
    gdb.printing.register_pretty_printer(gdb.current_objfile(), build_array_printer())
except ImportError:
    logger.warning("Could not register LvArray pretty printers.")


# import debugpy
# debugpy.listen(("0.0.0.0", 64018))
# debugpy.wait_for_client()
