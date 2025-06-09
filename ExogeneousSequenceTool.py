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

        # Add track_dim_reduction subcommand
        parser_track_dim_reduction = subparsers.add_parser("track_dim_reduction")
        SignalTrack.set_parser_track_dim_reduction(parser_track_dim_reduction)

        parser_mutagenesis = subparsers.add_parser("mutagenesis")
        Mutagenesis.set_parser_mutagenesis(parser_mutagenesis)

    @staticmethod
    def main(args: argparse.Namespace):
        if args.subcommand == "add_adapter":
            ExogeneousSequenceAssemble.add_adapter_main(args)
        elif args.subcommand == "concat":
            ExogeneousSequenceAssemble.concat_main(args)
        elif args.subcommand == "track_dim_reduction":
            SignalTrack.track_dim_reduction_main(args)
        elif args.subcommand == "mutagenesis":
            Mutagenesis.mutagenesis_main(args)
        else:
            raise ValueError(f"Subcommand {args.subcommand} not found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ExogeneousSequenceTool.set_parser(parser)
    args = parser.parse_args()

    ExogeneousSequenceTool.main(args)
