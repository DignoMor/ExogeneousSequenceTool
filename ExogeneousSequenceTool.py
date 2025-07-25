#!/usr/bin/env python

import argparse

from RGTools.ExogeneousSequences import ExogeneousSequences

from ExogeneousSequenceAssemble import ExogeneousSequenceAssemble
from SignalTrack import SignalTrack
from Mutagenesis import Mutagenesis
from Motif import Motif

class ExogeneousSequenceTool:
    @staticmethod
    def set_parser(parser: argparse.ArgumentParser):
        subparsers = parser.add_subparsers(dest="subcommand")
        
        parser_assemble = subparsers.add_parser("assemble", 
                                               help="Assemble exogeneous sequences.", 
                                               )
        ExogeneousSequenceAssemble.set_parser(parser_assemble)

        parser_track_dim_reduction = subparsers.add_parser("track_dim_reduction", 
                                                          help="Reduce the dimension of the signal track.", 
                                                          )
        SignalTrack.set_parser_track_dim_reduction(parser_track_dim_reduction)

        parser_mutagenesis = subparsers.add_parser("mutagenesis", 
                                                   help="Mutate the exogeneous sequences on a given base location.", 
                                                   )
        Mutagenesis.set_parser_mutagenesis(parser_mutagenesis)

        parser_gen_track = subparsers.add_parser("gen_track", 
                                                 help="Generate a signal track.", 
                                                 )
        SignalTrack.set_parser_gen_track(parser_gen_track)

        parser_motif_search = subparsers.add_parser("motif_search", 
                                                   help="Search for motifs in the exogeneous sequences.", 
                                                   )
        Motif.set_parser_motif_search(parser_motif_search)

    @staticmethod
    def main(args: argparse.Namespace):
        if args.subcommand == "assemble":
            ExogeneousSequenceAssemble.main(args)
        elif args.subcommand == "track_dim_reduction":
            SignalTrack.track_dim_reduction_main(args)
        elif args.subcommand == "mutagenesis":
            Mutagenesis.mutagenesis_main(args)
        elif args.subcommand == "gen_track":
            SignalTrack.gen_track_main(args)
        elif args.subcommand == "motif_search":
            Motif.motif_search_main(args)
        else:
            raise ValueError(f"Subcommand {args.subcommand} not found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ExogeneousSequenceTool.set_parser(parser)
    args = parser.parse_args()

    ExogeneousSequenceTool.main(args)
