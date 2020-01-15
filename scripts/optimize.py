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

    parser.add_argument("sequence", help='''DNA sequence to be optimized. Takes
        a single sequence as a string, or the name of a file with one or more 
        sequences in fasta or genbank format. If a file is provided, the 
        --intype argument has to be given as well.''')

    parser.add_argument("-i", "--intype", default=None, 
        help='''If `sequence` is a file, indicate the file type as `fasta` or 
        `genbank`.''')

    parser.add_argument("-o", "--outtype", default="fasta",
        help='''Save the resulting optimized sequences to a file, in the 
        indicated format. It can be `fasta`, `genbank` or `fasta/genbank` for 
        both.''')

    parser.add_argument("-v", "--verbose", help='''Show the constraints 
        evaluations, and optimization objectives score.''', action="count",
        default=0)

    parser.add_argument("-p", "--separate", help='''Save the optimized sequences
         from the input file in separate files, as opposed to all in one single 
        file''', action="store_true")

    parser.add_argument("-d", "--destination", default=os.getcwd(), 
        help='''Path for saving the resulting sequences. --outtype has to be 
        provided.''')

    parser.add_argument("-q", "--seqid", default="seq", 
        help='''Name/identifier of the sequence. Only used when a single
        sequence is provided as a string, and you want to save it to a file.''')

    # parser.add_argument("-m", "--modebt", default="optimize",
    #     help='''Mode to run back translation of the sequence, it can be 
    #     "optimize" to take into account codon frequencies, or "random" to give
    #     random codons (default=optimize).''')

    parser.add_argument("-j", "--jobname", default=None,
        help='''Name of the job that will go in the name of the results file. 
        --outtype has to be provided.''')

    parser.add_argument("-n", "--number", default=10, type=int, 
        help='''Number of optimized sequences to produce per each
        back-translated one (default=10).''')

    parser.add_argument("-b", "--nback", default=10, type=int, 
        help='''Number of random back-transalted sequences to be produced 
        (default=10). 0 to not back translate and optimize only from the 
        original sequence.''')

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


def get_jobname(jobname, sequence):
    """
    Gets the name for the output file

    Input
    -----
    jobname : str
        --jobname argument provided when calling the script, if not provided it
        is None
    sequence : str
        --sequence argument provided when calling the script
    """
    if jobname:
        return jobname
    elif os.path.isfile(sequence):
        # File name without the extension
        return os.path.splitext(os.path.basename(sequence))[0]
    else:
        return "coolercodonopt"


def save_file(optimized, destination, outtype, separate, jobname):
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
        jobname = get_jobname(args.jobname, args.sequence)
        
        if args.outtype == 'fasta/genbank':
            save_file(optimized, args.destination, 'fasta', args.separate,
                jobname)
            save_file(optimized, args.destination, 'genbank', args.separate,
                jobname)
        else:
            save_file(optimized, args.destination, args.outtype, args.separate,
                jobname)
