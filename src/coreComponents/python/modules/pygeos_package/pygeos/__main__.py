"""Command line tools for pygeos"""

import argparse
from pygeos import xml_processor


def main():
  """Entry point for the pygeos console script

  @arg input Input file name
  @arg -o/--output Output filename (default = randomly generated string)
  @arg -s/--schema GEOSX schema to use for validating the generated xml
  @arg -v/--verbose Verbosity level
  """

  # Parse the user arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('input', type=str, help='Input file name')
  parser.add_argument('-o', '--output', type=str, help='Output file name (otherwise, it is randomly genrated)', default='')
  parser.add_argument('-s', '--schema', type=str, help='GEOSX schema to use for validation', default='')
  parser.add_argument('-v', '--verbose', type=int, help='Verbosity of outputs', default=0)
  args = parser.parse_args()

  # Process the xml file
  output_name = xml_processor.process(args.input,
                                      outputFile=args.output,
                                      schema=args.schema,
                                      verbose=args.verbose)

  print(output_name)


if __name__ == "__main__":
  main()


