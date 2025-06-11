
import unittest
import argparse
import shutil
import os

from RGTools.ExogeneousSequences import ExogeneousSequences

from ExogeneousSequenceAssemble import ExogeneousSequenceAssemble

class ExogeneousSequenceAssembleTest(unittest.TestCase):
    def setUp(self):
        self.test_dir = "test_dir"
        os.makedirs(self.test_dir, exist_ok=True)

    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_add_adapter(self):
        args = argparse.Namespace(
            operation="add_adapter",
            fasta=os.path.join(self.test_dir, "test.fasta"),
            left_adapter_fasta=os.path.join(self.test_dir, "left_adapter.fasta"),
            right_adapter_fasta=os.path.join(self.test_dir, "right_adapter.fasta"),
            output_fasta=os.path.join(self.test_dir, "output.fasta"),
        )

        ExogeneousSequences.write_sequences_to_fasta(
            ["seq1", "seq2", "seq3"],
            ["ATCG", "TTGA", "CCAT"],
            args.fasta,
        )

        ExogeneousSequences.write_sequences_to_fasta(
            ["left_adapter"],
            ["AAA"],
            args.left_adapter_fasta,
        )

        ExogeneousSequences.write_sequences_to_fasta(
            ["right_adapter"],
            ["TTT"],
            args.right_adapter_fasta,
        )

        ExogeneousSequenceAssemble.main(args)

        output_es = ExogeneousSequences(args.output_fasta)
        self.assertEqual(output_es.get_region_bed_table().get_chrom_names()[1], "seq2")
        self.assertEqual(output_es.get_all_region_seqs()[1], "AAATTGATTT")

    def test_concat(self):
        args = argparse.Namespace(
            operation="concat",
            fasta5=os.path.join(self.test_dir, "fasta5.fasta"),
            fasta3=os.path.join(self.test_dir, "fasta3.fasta"),
            output_fasta=os.path.join(self.test_dir, "output.fasta"),
            id_method="5_3",
        )

        ExogeneousSequences.write_sequences_to_fasta(
            ["seq51", "seq52", "seq53"],
            ["ATCG", "TTGA", "CCAT"],
            args.fasta5,
        )
        
        ExogeneousSequences.write_sequences_to_fasta(
            ["seq31", "seq32", "seq33"],
            ["ATCG", "TTGA", "CCAT"],
            args.fasta3,
        )

        ExogeneousSequenceAssemble.main(args)

        output_es = ExogeneousSequences(args.output_fasta)
        self.assertEqual(output_es.get_region_bed_table().get_chrom_names()[1], "seq52_seq32")
        self.assertEqual(output_es.get_all_region_seqs()[1], "TTGATTGA")
        
        
if __name__ == "__main__":
    unittest.main()
