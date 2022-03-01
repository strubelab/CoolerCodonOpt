#!/usr/bin/env python3

"""
Script to parse the command line arguments and perform codon optimization
"""

from pathlib import Path
from datetime import datetime

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from parser_utils import parsing
from logger_utils import get_logger
from sequence_utils import save_file, read_sequences

from Optimizer import Optimizer


if __name__ == '__main__':

    logger = get_logger(__name__)

    try:
        args = parsing()
        sequences_dict = read_sequences(args.input)
    except ValueError:
        logger.error(
            'Please provide a valid file or directory name with the sequences in FASTA format.')
        raise

    for file, sequences in sequences_dict.items():
        logger.info(f'{file} : {len(sequences)} sequence(s) found.')
        optimizers = [Optimizer(seq) for seq in sequences]
        for opt in optimizers:
            opt.optimize()

            if args.verbose:
                if opt.is_nucleotide:
                    logger.info('ORIGINAL SEQUENCE:')
                    logger.info(f'> {opt.input_SeqRec.id}')
                    logger.info(f'{opt.input_SeqRec.seq}')
                    logger.info(opt.input_SeqRec.annotations["Constraints text summary"])
                    logger.info(opt.input_SeqRec.annotations["Objectives text summary"])

                logger.info('OPTIMIZED SEQUENCE:')
                logger.info(f'> {opt.best_optimized_SeqRec.id}')
                logger.info(f'{opt.best_optimized_SeqRec.seq}')
                logger.info(opt.best_optimized_SeqRec.annotations["Constraints text summary"])
                logger.info(opt.best_optimized_SeqRec.annotations["Objectives text summary"])

        try:
            save_file(file, optimizers, args.destination)
        except ValueError:
            logger.error('Please provide a valid name for the destination directory')
