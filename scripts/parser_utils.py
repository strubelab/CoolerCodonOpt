import argparse
from sequence_utils import validate_file

def parsing(args: list=None) -> argparse.Namespace:
    """
    Creates the argument parser instance and applies it to the command line
    input

    Input
    -----
    args : list
        List of the arguments to be parsed (only to be used for testing). If
        none is provided, it is taken from sys.argv
    """

    parser = argparse.ArgumentParser(description='''Takes a DNA sequence and
        optimizes it for expression in E. coli'''.replace('\n', ''))

    parser.add_argument("input", help='''Fasta file with the sequence(s)
        to be optimized, or a directory with fasta files.'''.replace('\n', ''),
        type=validate_file)
    
    parser.add_argument("-v", "--verbose", help='''Show the constraints 
        evaluations, and optimization objectives score.'''.replace('\n', ''),
        action="store_true")

    parser.add_argument("-d", "--destination", default=None, 
        help='''Path for saving the resulting sequences. It defaults to the same 
        directory as the input.'''.replace('\n', ''))

    return parser.parse_args(args)
