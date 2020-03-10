"""DFN Generation Scripts"""

import argparse
from pygeos.dfn import dfn_generator


def main():
  """Entry point for the pygeos console script

  @arg -i/--input Input file name
  @arg -o/--output Output filename
  @arg -n/--nIter Max number of iterations for DFN construction
  @arg -m/--margin Margin for DFN sets
  @arg -p/--percolation Calculate percolation for generated DFN
  """

  # Parse the user arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', type=str, help='Input file name')
  parser.add_argument('-o', '--output', type=str, help='Output file name')
  parser.add_argument('-n', '--nIter', type=int, help='Max number of iterations for DFN construction', default=10000)
  parser.add_argument('-m', '--margin', type=int, help='Margin for DFN sets', default=0.01)
  parser.add_argument('-p', '--percolation', type=int, help='Calculate percolation for generated DFN', default=True)
  args = parser.parse_args()

  dfn_generator.generate_from_xml(args.input,
                                  args.output,
                                  args.nIter,
                                  args.margin,
                                  args.percolation)


if __name__ == "__main__":
  main()


