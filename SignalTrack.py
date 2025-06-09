
import numpy as np

from RGTools.ExogeneousSequences import ExogeneousSequences

class SignalTrack:
    @staticmethod
    def set_parser_argmax(parser):
        ExogeneousSequences.set_parser_exogeneous_sequences(parser)

        parser.add_argument("--input_npy",
                            help="Path to the signal track npy file.",
                            required=True,
                            )

        parser.add_argument("--output_npy",
                            help="Path to the output stat npy file.",
                            required=True,
                            )
        
        parser.add_argument("--search_range",
                            help="Search range for the signal track. Format: 'start,end'. "
                                 "If not provided, the whole signal track will be used. "
                                 "Range follows half-open and 0-index convention.",
                            default=None, 
                            type=str,
                            )

    @staticmethod
    def argmax_main(args):
        es = ExogeneousSequences(args.fasta)
        es.load_region_anno_from_npy("signal_track", args.input_npy)
        signal_track = es.get_anno_arr("signal_track")

        if args.search_range: 
            start, end = map(int, args.search_range.split(","))
            signal_track[:, :start] = -np.inf
            signal_track[:, end:] = -np.inf

        output_stat = signal_track.argmax(axis=1, keepdims=True)
        es.load_region_anno_from_arr("stat", output_stat)
        es.save_anno_npy("stat", args.output_npy)
