from geosx_xml_tools.main import preprocess_parallel


def apply_xml_preprocessor():
    """
    Applies the xml preprocessor to the input file
    before handing it to GEOSX, and modifies the input
    arguments to point to the new file
    """
    return preprocess_parallel()
