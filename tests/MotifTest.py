
import os
import random
import shutil
import unittest
import argparse
import requests

import numpy as np

from RGTools.ExogeneousSequences import ExogeneousSequences
from RGTools.MemeMotif import MemeMotif

from Motif import Motif

class MotifTest(unittest.TestCase):
    def setUp(self) -> None:
        random.seed(76)
        
        self._test_path = "motif_test"
        self._meme_motif_path = os.path.join(self._test_path, "test_motifs.meme")
        self._fasta_path = os.path.join(self._test_path, "test_fasta.fasta")

        if not os.path.exists(self._test_path):
            os.makedirs(self._test_path)

        url = "https://meme-suite.org/meme/doc/examples/sample-dna-motif.meme"
        response = requests.get(url)

        # Save the file locally
        with open(self._meme_motif_path, "wb") as f:
            f.write(response.content)

        meme_motif = MemeMotif(self._meme_motif_path)
        pwm = meme_motif.get_motif_pwm("crp")
        alphabet = meme_motif.get_alphabet()

        seqs = []
        for _ in range(5):
            seq = "ATCG"
            # For each position in the PWM
            for pos in range(pwm.shape[0]):
                # Get probabilities for this position
                probs = pwm[pos]
                # Sample a nucleotide based on the probabilities
                nucleotide = random.choices(alphabet, weights=probs, k=1)[0]
                seq += nucleotide
            
            seq += "ATCG"
            seqs.append(seq)
        
        seqs[3] = seqs[3][:-3]

        if os.path.exists(self._fasta_path):
            os.remove(self._fasta_path)

        ExogeneousSequences.write_sequences_to_fasta(
            ["seq1", "seq2", "seq3", "seq4", "seq5"],
            seqs,
            self._fasta_path,
        )

        return super().setUp()
    
    def tearDown(self) -> None:
        if os.path.exists(self._test_path):
            shutil.rmtree(self._test_path)

        return super().tearDown()

    def test_motif_search(self):
        args = argparse.Namespace(
            fasta=self._fasta_path,
            motif_file=self._meme_motif_path,
            output_header=os.path.join(self._test_path, "test_motif_search"),
            estimate_background_freq=True,
            reverse_complement=False,
        )

        Motif.motif_search_main(args)
        crp_out = np.load(args.output_header + ".crp.npy")
        self.assertTrue((crp_out[:, 4] > 0).all())
        
    
