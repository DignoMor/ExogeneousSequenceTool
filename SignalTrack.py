
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

    @staticmethod
    def argmax_main(args):
        es = ExogeneousSequences(args.fasta)
        es.load_region_anno_from_npy("signal_track", args.input_npy)
        signal_track = es.get_anno_arr("signal_track")
        output_stat = signal_track.argmax(axis=1, keepdims=True)
        es.load_region_anno_from_arr("stat", output_stat)
        es.save_anno_npy("stat", args.output_npy)
