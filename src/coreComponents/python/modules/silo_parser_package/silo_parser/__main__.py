
from silo_parser.parse import parse_multiple
import argparse


def main():
  parser = argparse.ArgumentParser(description='Simple processor for GEOSX Silo files')

  parser.add_argument('-p',
                      '--plot_header',
                      type=str,
                      default='plot')

  parser.add_argument('-t',
                      '--type',
                      type=str,
                      default='Fracture_ElementFields')

  parser.add_argument('-v',
                      '--variables',
                      type=str,
                      nargs='+',
                      default=['elementAperture',
                               'effectiveAperture',
                               'elementArea',
                               'pressure',
                               'ruptureTime'])

  parser.add_argument('-n',
                      '--npar_proc',
                      type=int,
                      default=1)

  parser.add_argument('-f',
                      '--folder',
                      type=str,
                      default='data')

  args = parser.parse_args()

  parse_multiple(args.plot_header,
                 args.type,
                 args.variables,
                 parallel_folder=args.folder,
                 npar_proc=args.npar_proc)


if __name__ == "__main__":
  main()

