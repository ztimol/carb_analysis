import sys, os
import argparse
import pdb_reader
from config import get_analysis_config

parser = argparse.ArgumentParser(description="Calculate NOE from an md simulation")

parser.add_argument("-f", help="Select a config file", type=str)

args = parser.parse_args()


def main():
    env = {}
    env = {"results": {}}

    env["input_params"] = get_analysis_config(args.f)
    print(env)


if __name__ == "__main__":
    main()
