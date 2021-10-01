"""Command line tools for geosx_xml_tools"""

import argparse
from geosx_xml_tools import xml_processor


def main():
    """Entry point for the geosx_xml_tools console script

    @arg input Input file name
    @arg -o/--output Output filename (default = randomly generated string)
    @arg -s/--schema GEOSX schema to use for validating the generated xml
    @arg -v/--verbose Verbosity level
    """

    # Parse the user arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, action='append', help='Input file name')
    parser.add_argument('-o', '--output', type=str, help='Output file name (otherwise, it is randomly genrated)', default='')
    parser.add_argument('-s', '--schema', type=str, help='GEOSX schema to use for validation', default='')
    parser.add_argument('-v', '--verbose', type=int, help='Verbosity of outputs', default=0)
    parser.add_argument('-p', '--parameters', nargs='+', action='append', help='Parameter overrides', default=[])
    args = parser.parse_args()

    # Process the xml file
    output_name = xml_processor.process(args.input,
                                        outputFile=args.output,
                                        schema=args.schema,
                                        verbose=args.verbose,
                                        parameter_override=args.parameters)
    print(output_name)


if __name__ == "__main__":
    main()


