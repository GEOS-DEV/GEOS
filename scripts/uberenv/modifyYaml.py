import argparse
import os
import sys
import yaml

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument( "-i", "--input", type=os.PathLike, help = "input yaml file")
    parser.add_argument( "-o", "--output", type=os.PathLike, required = True, help="output yaml file" )
    parser.add_argument( "-y", "--yaml", type=str, action='append', required = True, help="list of (potentially nested) key-value pairs to add/overwrite in the yaml file." )
    args, _ = parser.parse_known_args()

    content = {}
    if args.input is not None:
      with open( args.input, 'r' ) as fin:
        content = yaml.load( fin, Loader=yaml.FullLoader )

    for ystr in args.yaml:
      ydict = yaml.safe_load( ystr )
      content = { **content, **ydict }

    with open( args.output, 'w' ) as fout:
       yaml.dump( content, fout )

    return 0

if __name__ == "__main__" and not sys.flags.interactive:
    sys.exit(main())
