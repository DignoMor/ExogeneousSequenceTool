#!/usr/bin/env python

import argparse

from RGTools.ExogeneousSequences import ExogeneousSequences

from ExogeneousSequenceAssemble import ExogeneousSequenceAssemble
from SignalTrack import SignalTrack
from Mutagenesis import Mutagenesis

class ExogeneousSequenceTool:
    @staticmethod
    def set_parser(parser: argparse.ArgumentParser):
        subparsers = parser.add_subparsers(dest="subcommand")
        parser_add_adapter = subparsers.add_parser("add_adapter")
        ExogeneousSequenceAssemble.set_parser_add_adapter(parser_add_adapter)

        parser_concat = subparsers.add_parser("concat")
        ExogeneousSequenceAssemble.set_parser_concat(parser_concat)

        parser_argmax = subparsers.add_parser("argmax")
        SignalTrack.set_parser_argmax(parser_argmax)

        parser_mutagenesis = subparsers.add_parser("mutagenesis")
        Mutagenesis.set_parser_mutagenesis(parser_mutagenesis)

    @staticmethod
    def main(args: argparse.Namespace):
        if args.subcommand == "add_adapter":
            ExogeneousSequenceAssemble.add_adapter_main(args)
        elif args.subcommand == "concat":
            ExogeneousSequenceAssemble.concat_main(args)
        elif args.subcommand == "argmax":
            SignalTrack.argmax_main(args)
        elif args.subcommand == "mutagenesis":
            Mutagenesis.mutagenesis_main(args)
        else:
            raise ValueError(f"Subcommand {args.subcommand} not found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ExogeneousSequenceTool.set_parser(parser)
    args = parser.parse_args()

    ExogeneousSequenceTool.main(args)
