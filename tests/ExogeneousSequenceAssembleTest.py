
import unittest
import argparse
import shutil
import os
import pandas as pd

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

    def test_barcode(self):
        # Create test files
        barcode_fasta = os.path.join(self.test_dir, "barcode.fasta")
        input_fasta1 = os.path.join(self.test_dir, "input1.fasta")
        input_fasta2 = os.path.join(self.test_dir, "input2.fasta")
        output_fasta = os.path.join(self.test_dir, "output.fasta")
        metadata_path = os.path.join(self.test_dir, "metadata.csv")

        # Write test sequences
        ExogeneousSequences.write_sequences_to_fasta(
            ["barcode1", "barcode2", "barcode3", "barcode4"],
            ["A", "T", "C", "G"],
            barcode_fasta,
        )

        ExogeneousSequences.write_sequences_to_fasta(
            ["seq1", "seq2"],
            ["ATCG", "TTGA"],
            input_fasta1,
        )

        ExogeneousSequences.write_sequences_to_fasta(
            ["seq3", "seq4"],
            ["CCAT", "GGTA"],
            input_fasta2,
        )

        args = argparse.Namespace(
            operation="barcode",
            barcode_fasta=barcode_fasta,
            input_fasta=[input_fasta1, input_fasta2],
            input_class=["class1", "class2"],
            output_fasta=output_fasta,
            metadata_path=metadata_path,
            barcode_method="5_3",
        )

        ExogeneousSequenceAssemble.main(args)

        # Check output fasta
        output_es = ExogeneousSequences(output_fasta)
        self.assertEqual(output_es.get_region_bed_table().get_chrom_names()[1], "seq2")
        self.assertEqual(output_es.get_all_region_seqs()[1], "TTTGAT")

        # Check metadata
        metadata_df = pd.read_csv(metadata_path)
        self.assertEqual(len(metadata_df), 4)  # 4 sequences total
        self.assertEqual(metadata_df.iloc[1]["barcode"], "T")
        self.assertEqual(metadata_df.iloc[1]["class"], "class1")
        self.assertEqual(metadata_df.iloc[1]["elem_id"], "seq2")
        self.assertEqual(metadata_df.iloc[1]["elem_seq"], "TTTGAT")

    def test_barcode_too_many_elements(self):
        # Create test files with more elements than barcodes
        barcode_fasta = os.path.join(self.test_dir, "barcode.fasta")
        input_fasta1 = os.path.join(self.test_dir, "input1.fasta")
        input_fasta2 = os.path.join(self.test_dir, "input2.fasta")
        output_fasta = os.path.join(self.test_dir, "output.fasta")
        metadata_path = os.path.join(self.test_dir, "metadata.csv")

        # Write test sequences
        ExogeneousSequences.write_sequences_to_fasta(
            ["barcode1", "barcode2"],  # Only 2 barcodes
            ["A", "T"],
            barcode_fasta,
        )

        ExogeneousSequences.write_sequences_to_fasta(
            ["seq1", "seq2"],
            ["ATCG", "TTGA"],
            input_fasta1,
        )

        ExogeneousSequences.write_sequences_to_fasta(
            ["seq3", "seq4"],
            ["CCAT", "GGTA"],
            input_fasta2,
        )

        args = argparse.Namespace(
            operation="barcode",
            barcode_fasta=barcode_fasta,
            input_fasta=[input_fasta1, input_fasta2],
            input_class=["class1", "class2"],
            output_fasta=output_fasta,
            metadata_path=metadata_path,
            barcode_method="5_3",
        )

        with self.assertRaises(ValueError) as context:
            ExogeneousSequenceAssemble.main(args)
        
        self.assertEqual(str(context.exception), 
                        "Total number of elements in the input fasta files is greater than the number of regions in the barcode fasta file.")
        
if __name__ == "__main__":
    unittest.main()
