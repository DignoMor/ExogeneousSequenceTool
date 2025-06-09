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
    
    def _create_test_data(self, args):
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
        return tracks

    def _get_test_args(self, search_range=None):
        return argparse.Namespace(
            fasta=os.path.join(self.test_dir, "fasta.fasta"),
            input_npy=os.path.join(self.test_dir, "input.npy"),
            output_npy=os.path.join(self.test_dir, "output.npy"),
            search_range=search_range,
        )

    def test_track_dim_reduction_argmax(self):
        args = self._get_test_args()
        args.operation = "argmax"
        tracks = self._create_test_data(args)
        
        SignalTrack.track_dim_reduction_main(args)

        output_track = np.load(args.output_npy)
        self.assertEqual(output_track.shape, (3, 1))
        self.assertEqual(output_track[0, 0], 7)  # max at index 7
        self.assertEqual(output_track[1, 0], 4)  # max at index 4
        self.assertEqual(output_track[2, 0], 5)  # max at index 5

    def test_track_dim_reduction_max(self):
        args = self._get_test_args()
        args.operation = "max"
        tracks = self._create_test_data(args)
        
        SignalTrack.track_dim_reduction_main(args)

        output_track = np.load(args.output_npy)
        self.assertEqual(output_track.shape, (3, 1))
        self.assertEqual(output_track[0, 0], 3)  # max value is 3
        self.assertEqual(output_track[1, 0], 1)  # max value is 1
        self.assertEqual(output_track[2, 0], 2)  # max value is 2

    def test_track_dim_reduction_argmin(self):
        args = self._get_test_args()
        args.operation = "argmin"
        tracks = self._create_test_data(args)
        
        SignalTrack.track_dim_reduction_main(args)

        output_track = np.load(args.output_npy)
        self.assertEqual(output_track.shape, (3, 1))
        self.assertEqual(output_track[0, 0], 0)  # min at index 0
        self.assertEqual(output_track[1, 0], 0)  # min at index 0
        self.assertEqual(output_track[2, 0], 0)  # min at index 0

    def test_track_dim_reduction_min(self):
        args = self._get_test_args()
        args.operation = "min"
        tracks = self._create_test_data(args)
        
        SignalTrack.track_dim_reduction_main(args)

        output_track = np.load(args.output_npy)
        self.assertEqual(output_track.shape, (3, 1))
        self.assertEqual(output_track[0, 0], 0)  # min value is 0
        self.assertEqual(output_track[1, 0], 0)  # min value is 0
        self.assertEqual(output_track[2, 0], 0)  # min value is 0

    def test_track_dim_reduction_with_search_range(self):
        args = self._get_test_args(search_range="2,6")
        args.operation = "argmax"
        tracks = self._create_test_data(args)
        
        SignalTrack.track_dim_reduction_main(args)

        output_track = np.load(args.output_npy)
        self.assertEqual(output_track.shape, (3, 1))
        self.assertEqual(output_track[0, 0], 5)  # max in range at index 5
        self.assertEqual(output_track[1, 0], 4)  # max in range at index 4
        self.assertEqual(output_track[2, 0], 5)  # max in range at index 5




