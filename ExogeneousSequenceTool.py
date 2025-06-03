#!/usr/bin/env python

import argparse

from RGTools.ExogeneousSequences import ExogeneousSequences

from ExogeneousSequenceAssemble import ExogeneousSequenceAssemble

class ExogeneousSequenceTool:
    @staticmethod
    def set_parser(parser: argparse.ArgumentParser):
        subparsers = parser.add_subparsers(dest="subcommand")
        parser_add_adapter = subparsers.add_parser("add_adapter")
        ExogeneousSequenceAssemble.set_parser_add_adapter(parser_add_adapter)

    @staticmethod
    def main(args: argparse.Namespace):
        if args.subcommand == "add_adapter":
            ExogeneousSequenceAssemble.add_adapter_main(args)
        else:
            raise ValueError(f"Subcommand {args.subcommand} not found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ExogeneousSequenceTool.set_parser(parser)
    args = parser.parse_args()

    ExogeneousSequenceTool.main(args)
