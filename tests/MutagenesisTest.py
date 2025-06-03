
import unittest
import argparse
import shutil
import os

import numpy as np

from RGTools.ExogeneousSequences import ExogeneousSequences
from Mutagenesis import Mutagenesis

class MutagenesisTest(unittest.TestCase):
    def setUp(self):
        self.test_dir = "test_dir"
        os.makedirs(self.test_dir, exist_ok=True)
        super().setUp()
    
    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
        super().tearDown()
    
    def test_mutagenesis(self):
        args = argparse.Namespace(
            fasta=os.path.join(self.test_dir, "fasta.fasta"),
            loc_npy=os.path.join(self.test_dir, "loc.npy"),
            mut_fasta=os.path.join(self.test_dir, "mut.fasta"),
            output_fasta=os.path.join(self.test_dir, "output.fasta"),
        )

        ExogeneousSequences.write_sequences_to_fasta( 
            ["seq1", "seq2", "seq3"],
            ["ATCG", "TTGA", "CCAT"],
            args.fasta,
        )

        ExogeneousSequences.write_sequences_to_fasta(
            ["mut1", "mut2", "mut3"],
            ["C", "A", "T"],
            args.mut_fasta,
        )

        loc_arr = np.array([0, 1, 2]).reshape(-1, 1)
        np.save(args.loc_npy, loc_arr)

        Mutagenesis.mutagenesis_main(args)

        output_es = ExogeneousSequences(args.output_fasta)
        self.assertEqual(output_es.get_region_bed_table().get_chrom_names()[1], "seq2_mut_1_A")
        self.assertEqual(output_es.get_all_region_seqs()[1], "TAGA")


