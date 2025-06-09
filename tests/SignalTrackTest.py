
import unittest
import argparse
import shutil
import os

import numpy as np

from RGTools.ExogeneousSequences import ExogeneousSequences
from SignalTrack import SignalTrack

class SignalTrackTest(unittest.TestCase):
    def setUp(self):
        self.test_dir = "test_dir"
        os.makedirs(self.test_dir, exist_ok=True)
        super().setUp()

    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
        super().tearDown()
    
    def test_argmax(self):
        args = argparse.Namespace(
            fasta=os.path.join(self.test_dir, "fasta.fasta"),
            input_npy=os.path.join(self.test_dir, "input.npy"),
            output_npy=os.path.join(self.test_dir, "output.npy"),
            search_range=None,
        )
        
        ExogeneousSequences.write_sequences_to_fasta(
            ["seq1", "seq2", "seq3"],
            ["ATCG", "TTGA", "CCAT"],
            args.fasta,
        )

        tracks = np.zeros((3, 10))
        tracks[0, 3:7] = 1
        tracks[1, 8] = 1
        tracks[2, 4] = 1
        np.save(args.input_npy, tracks)
        
        SignalTrack.argmax_main(args)

        output_track = np.load(args.output_npy)
        self.assertEqual(output_track.shape, (3, 1))
        self.assertEqual(output_track[0, 0], 3)
        self.assertEqual(output_track[1, 0], 8)
        self.assertEqual(output_track[2, 0], 4)

    def test_argmax_with_search_range(self):
        args = argparse.Namespace(
            fasta=os.path.join(self.test_dir, "fasta.fasta"),
            input_npy=os.path.join(self.test_dir, "input.npy"),
            output_npy=os.path.join(self.test_dir, "output.npy"),
            search_range="2,6",
        )
        
        ExogeneousSequences.write_sequences_to_fasta(
            ["seq1", "seq2", "seq3"],
            ["ATCG", "TTGA", "CCAT"],
            args.fasta,
        )

        tracks = np.zeros((3, 10))
        tracks[0, 3:7] = 1
        tracks[0, 5] = 2
        tracks[0, 7] = 3
        tracks[1, 4] = 1
        tracks[2, 4] = 1
        tracks[2, 5] = 2
        np.save(args.input_npy, tracks)
        
        SignalTrack.argmax_main(args)

        output_track = np.load(args.output_npy)
        self.assertEqual(output_track.shape, (3, 1))
        self.assertEqual(output_track[0, 0], 5)
        self.assertEqual(output_track[1, 0], 4)
        self.assertEqual(output_track[2, 0], 5)




