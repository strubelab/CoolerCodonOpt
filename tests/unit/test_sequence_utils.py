import unittest
import os, sys, tempfile, shutil
from pathlib import Path, PosixPath
from Bio.SeqRecord import SeqRecord
import argparse
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, 'scripts'))

from sequence_utils import validate_file, read_sequences, check_nucleotide, \
    get_codon_probs, back_translate, detect_ATA_ATA, save_file

class TestValidateFile(unittest.TestCase):
    """
    Test the file validation function
    """

    def test_validate_right_input(self) -> None:
        input = os.path.join(os.path.dirname(__file__), 'fixtures/example.fasta')
        inp = validate_file(input)
        self.assertEqual(inp, Path(input))
    
    def test_validate_incorrect_input(self) -> None:
        input = 'false_input.fasta'
        with self.assertRaises(ValueError):
            inp = validate_file(input)

class TestReadSequences(unittest.TestCase):
    """
    Test the sequence reading function
    """
    def setUp(self) -> None:
        self.input_file = Path(os.path.join(os.path.dirname(__file__), 'fixtures/example.fasta'))
        self.input_dir = Path(os.path.join(os.path.dirname(__file__), 'fixtures'))
        self.input_file2 = Path(os.path.join(os.path.dirname(__file__), 'fixtures/example2.fasta'))


    def test_read_single_file(self) -> None:
        """
        A single fasta file is read
        """
        sequences = read_sequences(self.input_file)
        self.assertTrue(isinstance(sequences, dict))
        self.assertTrue(len(sequences)==1)
        
        result_file = list(sequences.keys())[0]
        self.assertTrue(isinstance(result_file, PosixPath))
        self.assertTrue(str(result_file)==str(self.input_file))

        read_SeqRecords = list(sequences.values())[0]
        self.assertTrue(len(read_SeqRecords)==2)
        self.assertTrue(isinstance(read_SeqRecords[0], SeqRecord))

    def test_read_directory(self) -> None:
        """
        A directory with fasta sequences is read
        """
        sequences = read_sequences(self.input_dir)
        self.assertTrue(isinstance(sequences, dict))
        self.assertTrue(len(sequences)==2)

        result_file0 = list(sequences.keys())[0]
        self.assertTrue(isinstance(result_file0, PosixPath))

        result_files = [str(p) for p in list(sequences.keys())]
        input_files = [str(p) for p in [self.input_file, self.input_file2]]
        self.assertTrue(result_files[0] in input_files)
        self.assertTrue(result_files[1] in input_files)

        read_SeqRecords = list(sequences.values())
        self.assertTrue(len(read_SeqRecords[0])==2)
        self.assertTrue(len(read_SeqRecords[1])==2)
        self.assertTrue(isinstance(read_SeqRecords[0][0], SeqRecord))


class TestCheckNucleotide(unittest.TestCase):
    """
    Test the check_nucleotide method
    """
    def test_nucleotide_sequence(self) -> None:
        """
        Give a nucleotide sequence as input
        """
        seq = 'ACTGACTGACTGACTGACTG'
        is_nucleotide = check_nucleotide(seq)
        self.assertTrue(is_nucleotide)

    def test_aminoacid_sequence(self) -> None:
        """
        Give an amino acid sequence as input
        """
        seq = 'ACDEFGHIKLMNPQ'
        is_nucleotide = check_nucleotide(seq)
        self.assertFalse(is_nucleotide)
    

class TestGetCodonProbs(unittest.TestCase):
    """
    Test the get_codon_probs method
    """
    def test_get_codons(self) -> None:
        """
        Test that get_codon_probs() returns a dictionary of the form
        {'*': [('TAA', 'TAG', 'TGA'), (0.64, 0.07, 0.29)], 'A': ... }
        """
        codon_probs = get_codon_probs()

        self.assertTrue(len(codon_probs)==21)
        self.assertTrue(len(codon_probs['A'])==2)
        self.assertTrue(len(codon_probs['A'][0])==4)
        self.assertTrue(len(codon_probs['A'][1])==4)
        self.assertTrue(isinstance(codon_probs['A'][0], np.ndarray))
        self.assertTrue(isinstance(codon_probs['A'][1], np.ndarray))

class TestBackTranslate(unittest.TestCase):
    """
    Test the back_translate method

    Args:
        unittest ([type]): [description]
    """
    def setUp(self) -> None:
        """
        Setup random number generator
        """
        entropy = 317753445747471334433479883993670272010
        sq = np.random.SeedSequence(entropy)
        self.rng = np.random.default_rng(sq)

        # E probabilities -> ['GAA', 'GAG'], [0.69, 0.31]
        # C probabilities -> ['TGC', 'TGT'], [0.56, 0.44])
        self.aa_sequence = 'E'*10 + 'C'*10
    
    def test_back_translate_optimized(self) -> None:
        """
        Test back translation with optimized frequencies
        """
        back_translated = back_translate(self.aa_sequence, rng=self.rng)
        # 7 GAA, 3 GAG, 5 TGC, 5 TGT
        test_back_translated = 'GAGGAAGAAGAGGAAGAAGAAGAGGAAGAATGCTGTTGTTGTTGCTGTTGCTGTTGCTGC'
        self.assertEqual(back_translated, test_back_translated)

    def test_back_translate_random(self) -> None:
        """
        Test back translation with random frequencies
        """
        back_translated = back_translate(self.aa_sequence, mode='random', rng=self.rng)
        # 4 GAA, 6 GAG, 5 TGC, 5 TGT
        test_back_translated = 'GAGGAGGAAGAAGAGGAAGAGGAGGAGGAATGCTGCTGTTGCTGTTGTTGTTGTTGCTGC'
        self.assertEqual(back_translated, test_back_translated)


class TestDetectATAATA(unittest.TestCase):
    """
    Test the detect_ATA_ATA method

    Args:
        unittest ([type]): [description]
    """

    def test_ATA_ATA(self):
        """
        Test case where two consecutive ATA codons are present
        """
        seq = 'GGGCCCATAATATTTGGG'
        self.assertTrue(detect_ATA_ATA(seq))
    
    def test_not_ATA_ATA(self):
        """
        No consecutive ATA codons
        """
        seq = 'GGGCCCAAATTTAAACCC'
        self.assertFalse(detect_ATA_ATA(seq))

