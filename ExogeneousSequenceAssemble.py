
from RGTools.ExogeneousSequences import ExogeneousSequences

class ExogeneousSequenceAssemble:
    @staticmethod
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

    @staticmethod
    def set_parser_concat(parser):
        parser.add_argument("--fasta5", 
                            help="Path to the 5' fasta file.",
                            required=True,
                            )

        parser.add_argument("--fasta3", 
                            help="Path to the 3' fasta file.",
                            required=True,
                            )

        parser.add_argument("--output_fasta", 
                            help="Path to the output fasta file.",
                            required=True,
                            default=None,
                            )
        
        parser.add_argument("--id_method",
                            help="Method to use to generate the id of the output fasta file.",
                            choices=["5", "3", "5_3"], 
                            default="5_3",
                            )

    @staticmethod
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

    @staticmethod
    def concat_main(args):
        fasta5_es = ExogeneousSequences(args.fasta5)
        fasta3_es = ExogeneousSequences(args.fasta3)

        output_seqs = []
        output_seq_ids = []
        for region5, region3 in zip(fasta5_es.get_region_bed_table().iter_regions(), fasta3_es.get_region_bed_table().iter_regions()):
            seq5 = fasta5_es.get_region_seq(region5["chrom"], 
                                            region5["start"], 
                                            region5["end"], 
                                            )
            seq3 = fasta3_es.get_region_seq(region3["chrom"], 
                                            region3["start"], 
                                            region3["end"], 
                                            )
            oseq = seq5 + seq3
            if args.id_method == "5":
                oid = region5["chrom"]
            elif args.id_method == "3":
                oid = region3["chrom"]
            elif args.id_method == "5_3":
                oid = region5["chrom"] + "_" + region3["chrom"]
            else:
                raise ValueError(f"Invalid id method: {args.id_method}")
            
            output_seqs.append(oseq)
            output_seq_ids.append(oid)

        ExogeneousSequences.write_sequences_to_fasta(output_seq_ids, output_seqs, args.output_fasta)

