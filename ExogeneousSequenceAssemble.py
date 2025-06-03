
from RGTools.ExogeneousSequences import ExogeneousSequences

class ExogeneousSequenceAssemble:
    def set_parser_add_adapter(parser):
        ExogeneousSequences.set_parser_exogeneous_sequences(parser)

        parser.add_argument("--left_adapter_fasta", 
                            help="Path to the 5' adapter fasta file.",
                            required=True,
                            )
        
        parser.add_argument("--right_adapter_fasta", 
                            help="Path to the 3' adapter fasta file.",
                            default=None,
                            )
        
        parser.add_argument("--output_fasta", 
                            help="Path to the output fasta file.",
                            required=True,
                            default=None,
                            )


    def add_adapter_main(args):
        
        if args.left_adapter_fasta:
            left_adapter_seqs = ExogeneousSequences(args.left_adapter_fasta).get_all_region_seqs()
            if len(left_adapter_seqs) != 1:
                raise ValueError("Left adapter fasta file must contain exactly one sequence.")
            left_adapter_seq = left_adapter_seqs[0]
        else:
            left_adapter_seq = ""

        if args.right_adapter_fasta:
            right_adapter_seqs = ExogeneousSequences(args.right_adapter_fasta).get_all_region_seqs()
            if len(right_adapter_seqs) != 1:
                raise ValueError("Right adapter fasta file must contain exactly one sequence.")
            right_adapter_seq = right_adapter_seqs[0]
        else:
            right_adapter_seq = ""

        input_es = ExogeneousSequences(args.fasta)

        output_seqs = []
        output_seq_ids = []
        for region in input_es.get_region_bed_table().iter_regions():
            region_seq = input_es.get_region_seq(region["chrom"], 
                                                 region["start"], 
                                                 region["end"], 
                                                 )
            output_seqs.append(left_adapter_seq + region_seq + right_adapter_seq)
            output_seq_ids.append(region["chrom"])

        ExogeneousSequences.write_sequences_to_fasta(output_seq_ids, output_seqs, args.output_fasta)
