import numpy as np

from RGTools.ExogeneousSequences import ExogeneousSequences

class SignalTrack:
    @staticmethod
    def _set_track_dim_reduction_common_parser_args(parser):
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
    def set_parser_track_dim_reduction(parser):
        """Set up parser for track dimension reduction operations."""
        subparsers = parser.add_subparsers(dest="operation", required=True)
        
        parser_max = subparsers.add_parser("max")
        SignalTrack._set_track_dim_reduction_common_parser_args(parser_max)
        
        parser_argmax = subparsers.add_parser("argmax")
        SignalTrack._set_track_dim_reduction_common_parser_args(parser_argmax)
        
        parser_min = subparsers.add_parser("min")
        SignalTrack._set_track_dim_reduction_common_parser_args(parser_min)
        
        parser_argmin = subparsers.add_parser("argmin")
        SignalTrack._set_track_dim_reduction_common_parser_args(parser_argmin)

    @staticmethod
    def _track_dim_reduction(args, operation):
        es = ExogeneousSequences(args.fasta)
        es.load_region_anno_from_npy("signal_track", args.input_npy)
        signal_track = es.get_anno_arr("signal_track")

        if args.search_range: 
            start, end = map(int, args.search_range.split(","))
            signal_track[:, :start] = -np.inf
            signal_track[:, end:] = -np.inf

        if operation == "argmax":
            output_stat = signal_track.argmax(axis=1, keepdims=True)
        elif operation == "max":
            output_stat = signal_track.max(axis=1, keepdims=True)
        elif operation == "argmin":
            output_stat = signal_track.argmin(axis=1, keepdims=True)
        elif operation == "min":
            output_stat = signal_track.min(axis=1, keepdims=True)
        else:
            raise ValueError(f"Unknown operation: {operation}")

        es.load_region_anno_from_arr("stat", output_stat)
        es.save_anno_npy("stat", args.output_npy)

    @staticmethod
    def track_dim_reduction_main(args):
        SignalTrack._track_dim_reduction(args, args.operation)
