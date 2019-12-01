"""
Script to parse the command line arguments and perform codon optimization
"""

import argparse, os
import codonopt


def parsing(args=None):
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
        optimizes it for expression in E. coli''')

    parser.add_argument("sequence", help='''DNA sequence to be optimized. Takes
        a single sequence as a string, or a text file with FASTA sequences.''')

    parser.add_argument("-v", "--verbose", help='''Show the constraints 
        evaluations, and optimization objectives score.''', action="count",
        default=0)

    parser.add_argument("-b", "--nback", default=10, type=int, 
        help='''Number of random back-transalted sequences to be produced 
        (default=10). 0 to not back translate and optimize only from the 
        original sequence.''')

    parser.add_argument("-n", "--number", default=10, type=int, 
        help='''Number of optimized sequences to produce per each
        back-translated one (default=10).''')

    parser.add_argument("-s", "--save", help='''Save the resulting optimized
        sequences to a file (only the highest scoring).''', action="store_true")

    parser.add_argument("-e", "--saveall", help='''Save the best optimized 
        sequences from each back-translation for comparison.''', 
        action="store_true")

    parser.add_argument("-d", "--destination", default=os.getcwd(), 
        help='''Directory path for saving the resulting sequences and scores, 
        if --save or --saveall options are selected.''')

    parser.add_argument("-a", "--atg", help='''Add ATG at the beginning of the
        sequence.''', action="store_true")

    parser.add_argument("-q", "--seqname", default="seq", 
        help='''Name/identifier of the sequence. Only used when a single
        sequence is provided.''')

    parser.add_argument("-m", "--modebt", default="optimize",
        help='''Mode to run back translation of the sequence, it can be 
        "optimize" to take into account codon frequencies, or "random" to give
        random codons (default=optimize).''')

    parser.add_argument("-j", "--jobname", default="coolercodonopt",
        help='''Name of the job that will go in the name of the results file 
        when activating --save.''')

    return parser.parse_args(args)


if __name__ == '__main__':

    args = parsing()
    sequence = 'ATG'+args.sequence if args.atg else args.sequence

    print('')
    print('Producing %d sequences for each input sequence...\n' % \
        (args.nback*args.number or args.number))
    
    codonopt.optimize(sequence, args.nback, args.number, args.saveall,
        args.destination, args.verbose, args.seqname, args.save, 
        args.modebt, args.jobname)