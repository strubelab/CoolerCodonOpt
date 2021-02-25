#!/usr/bin/env python3

"""
Script to parse the command line arguments and perform codon optimization
"""

import argparse, os
import codonopt

from Bio import SeqIO


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

    parser.add_argument("sequence", help='''DNA or protein sequence to be optimized. 
        Takes a single sequence as a string, the name of a file with one or more 
        sequences in fasta or genbank format, or a directory with the sequences 
        (the termination for fasta files must be `.fasta` or `.fa`, and `.gbk` for 
        genbank). If a file is provided, the --intype argument has to be given as well.
        '''.replace('\n', ''))

    parser.add_argument("-i", "--intype", default=None, 
        help='''If `sequence` is a file or a directory with files, indicate the 
        file type as `fasta` or `genbank`.'''.replace('\n', ''))

    parser.add_argument("-o", "--outtype", default="fasta",
        help='''Save the resulting optimized sequences to a file, in the 
        indicated format. It can be `fasta`, `genbank` or `fasta/genbank` for 
        both.'''.replace('\n', ''))

    parser.add_argument("-v", "--verbose", help='''Show the constraints 
        evaluations, and optimization objectives score.'''.replace('\n', ''),
        action="count", default=0)

    parser.add_argument("-p", "--separate", help='''Save the optimized sequences
         from the input file in separate files, as opposed to all in one single 
        file'''.replace('\n', ''), action="store_true")

    parser.add_argument("-d", "--destination", default=os.getcwd(), 
        help='''Path for saving the resulting sequences. --outtype has to be 
        provided.'''.replace('\n', ''))

    parser.add_argument("-q", "--seqid", default="seq", 
        help='''Name/identifier of the sequence. Only used when a single
        sequence is provided as a string, and you want to save it to a file.''')

    return parser.parse_args(args)


def print_verbose(og_sequence, opt_sequence):
    """
    Print the constraints evaluations and objectives score for the best sequence
    
    Input
    -----
    og_sequence : SeqRecord
        Original sequence
    opt_sequence : SeqRecord
        Optimized sequence
    """
    if og_sequence:
        print('> Original sequence')
        print(str(og_sequence.seq) + '\n')
        print(og_sequence.annotations["Constraints text summary"])
        print(og_sequence.annotations["Objectives text summary"])

    print('> Optimized sequence')
    print(str(opt_sequence.seq) + '\n')
    print(opt_sequence.annotations["Constraints text summary"])
    print(opt_sequence.annotations["Objectives text summary"])


def message(verbose, originals, optimized):
    """
    Prints the verbose messages to stdout if they are required

    Input
    -----
    verbose : int
        Count of the number of -v flags provided when calling the script
    originals : list
        List of SeqRecord objects with the original sequences if they are DNA,
        else the list contains None objects
    optimized : list
        List of SeqRecord objects with the optimzied sequences
    """
    if verbose==1:
        for i in range(len(optimized)):
            print_verbose(originals[i], optimized[i])
    elif verbose > 1:
        for i in range(len(optimized)):
            print('Optimized score: %.4f' % \
                optimized[i].annotations["Objectives score"])


def save_file(optimized, destination, outtype, separate, input_seq):
    """
    Save the optimized sequence(s) to a file

    Input
    -----
    optimized : list
        List of SeqRecord objects with the optimized sequences
    destination : str
        --destination argument
    outtype : str
        --outtype argument
    separate : boolean
        --separate argument, whether or not to save the optimzied sequences in
        individual files, default = False
    jobnale : str
        Name of the output file name
    """
    if os.path.isfile(input_seq):
        fname = os.path.splitext(os.path.basename(input_seq))[0]
    else:
        fname = "coolercodonopt"

    extension = "_optimized.fasta" if outtype == "fasta" else "_optimized.gbk"

    if separate:
        for opt in optimized:
            fname = os.path.join(destination, opt.name+extension)
            SeqIO.write(opt, fname, outtype)
    else:
        fname = os.path.join(destination, jobname+extension)
        SeqIO.write(optimized, fname, outtype)


if __name__ == '__main__':

    args = parsing()

    print('')
    print('Choosing the best out of %d optimized sequences for each input ' % \
        (args.nback*args.number or args.number) + 'sequence...\n')
    
    originals, optimized = codonopt.optimize(args.sequence, args.intype, 
        args.nback, args.number, args.seqid)

    if args.verbose:
        message(args.verbose, originals, optimized)

    if args.outtype:
        
        if args.outtype == 'fasta/genbank':
            save_file(optimized, args.destination, 'fasta', args.separate,
                args.sequence)
            save_file(optimized, args.destination, 'genbank', args.separate,
                args.sequence)
        else:
            save_file(optimized, args.destination, args.outtype, args.separate,
                jobname)
