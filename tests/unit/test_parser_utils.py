import unittest
import os, sys
from pathlib import Path
import argparse

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, 'scripts'))

from parser_utils import parsing

class TestParserUtils(unittest.TestCase):
    """
    Class for testing the parser utils
    """
    def setUp(self) -> None:
        """
        Setup the arguments to be tested
        """
        self.input = os.path.join(os.path.dirname(__file__), 'fixtures/example.fasta')
        self.dest = os.path.join(os.path.dirname(__file__), 'fixtures')
        self.argstring = [self.input, '--verbose', '--destination', self.dest]
        self.args = parsing(self.argstring)

    def test_input_path_parsed(self) -> None:
        """
        Test that the input is parsed correctly
        """
        self.assertEqual(self.args.input, Path(self.input))

    def test_verbose_parsed(self) -> None:
        """
        Test that the value of the verbose argument is captured correctly
        """
        self.assertEqual(self.args.verbose, True)

    def test_destination_parsed(self) -> None:
        """
        Test that the destination directory was parsed
        """
        self.assertEqual(self.args.destination, self.dest)
        
