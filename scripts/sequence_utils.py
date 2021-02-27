"""
Miscelaneous functions
"""
import numpy as np
from pathlib import Path

import python_codon_tables
from Bio.Data import IUPACData
from Bio import SeqIO
from numpy.random import Generator


def validate_file(input: str) -> bool:
    """
    Validate that input is an existing file or directory

    Args:
        input (str): input file or directory

    Returns:
        bool
    """
    inp = Path(input)
    if not inp.exists():
        raise ValueError

    return inp


def read_sequences(input: Path) -> dict:
    """
    Validate that input is an existing file or directory, and read the sequences

    Args:
        inp (str): input fasta file or directory
    
    Returns:
        dictionary with { fname1:[sequences], fname2:[sequences] ... } for all
        files in the input
    """
    if input.is_file():
        return {input : list(SeqIO.parse(input, 'fasta'))}
    elif input.is_dir():
        flist = [list(input.glob(pattern)) for pattern in ['*.fasta', '*.fa', '*.faa']]
        flist = [f for sublist in flist for f in sublist ]
        return {f : list(SeqIO.parse(f, 'fasta')) for f in flist} 


def check_nucleotide(sequence: str) -> bool:
    """
    Checks wether a given sequence is nucleotide

    Args:
        sequence (str)

    Returns:
        bool
    """
    amino_acids = (IUPACData.protein_letters
                    .replace('A','')
                    .replace('C','')
                    .replace('G','')
                    .replace('T',''))
    
    for aa in amino_acids:
        if aa in sequence:
            # logger.info(f'{sequence.id} : Protein sequence - {len(sequence)} amino acids.')
            return False
    
    # logger.info(f'{sequence.id} : DNA sequence - {len(sequence)} nucleotides.')
    return True


def get_codon_probs():
    """
    Get dictionary with codon probabilities, in the form:
    {'*': [('TAA', 'TAG', 'TGA'), (0.64, 0.07, 0.29)], 'A': ... }
    """
    # Dictionary of the form
    # {'*': {'TAA': 0.64, 'TAG': 0.07, 'TGA': 0.29}, 'A': ... }
    codon_table = python_codon_tables.get_codons_table("e_coli")

    # A modification is necessary since the probabilities of G do not
    # sum 1 ... substracted 0.01 to the value
    codon_table['G']['GGT'] = 0.33

    # Change it to
    # {'*': [('TAA', 'TAG', 'TGA'), (0.64, 0.07, 0.29)], 'A': ... }
    codons_probabilities = {}
    for aa, codons_ps in codon_table.items():
        codons, ps = zip(*codons_ps.items())
        codons_probabilities[aa] = [np.array(codons), np.array(ps)]

    return codons_probabilities


def back_translate(sequence:str, mode:str='optimize', rng:Generator=None) -> list:
    """
    Generate the back translations of the nucleotide sequences, either in
    'optimized' mode following the frequencies of a reference organism
    (default E. coli), or 'random' with equal frequencies for all the codons.

    Args:
        sequence (str): Amino acid sequence to back translate
        mode (str, optional): Back translation mode. Defaults to 'optimize'.
        rng (Generator, optional): random number generator, for testing purposes

    Returns:
        list: List of strings of back translated sequences.
    """
    # Initialize numpy's default random number generator
    rng = rng or np.random.default_rng()
    
    codons_probabilities = get_codon_probs()
    bt_seq = ''

    if mode == 'optimize':
        for aa in sequence:
            codons, ps = codons_probabilities[aa]
            bt_seq += rng.choice(codons, p=ps)
    elif mode == 'random':
        for aa in sequence:
            codons = codons_probabilities[aa][0]
            bt_seq += rng.choice(codons)

    return bt_seq


def detect_ATA_ATA(sequence:str) -> bool:
    """
    Detects the presence of two consecutive ATA codons in the given DNA sequence

    Input
    -----
    sequence : str
    """

    codon_positions = list(range(0,len(sequence), 3))
    
    for index in codon_positions:
        icodon = sequence[index:index+3]
        
        if icodon == 'ATA':

            if index+3 < len(sequence):
                jcodon = sequence[index+3:index+6]
            
                if jcodon == 'ATA':
                    return True

    return False


def save_file(file:Path, optimizers:list, destination:str=None):
    """
    Writes the optimized sequences to a file

    Args:
        file (Path): Input file.
        optimizers (list): List of Optimizer objects with the optimized sequences.
        destination (str): String with the destination directory name.

    Raises:
        ValueError: Raised if the destination directory is not valid.
    """
    if destination:
        dest = Path(destination)
        if not dest.is_dir():
            raise ValueError

        f_out = dest / (file.stem + '_optimized.fasta')
    else:
        f_out = file.parent / (file.stem + '_optimized.fasta')
    
    optimized_SeqRecs = [opt.best_optimized_SeqRec for opt in optimizers]
    
    SeqIO.write(optimized_SeqRecs, f_out, 'fasta')