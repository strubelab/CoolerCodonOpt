"""
Script to perform codon optimization using the dnachisel library
"""

import os, random, copy

import dnachisel as dn
import dnachisel.builtin_specifications as spec
import python_codon_tables as pct
import numpy as np

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

from operator import itemgetter


# def get_sequences(sequence):
#     """
#     Extracts the sequences from a fasta or text file and saves them to a list

#     Input
#     -----
#     sequence : str
#         Name of fasta or text file

#     Output
#     -----
#     sequences: list
#         List of tuples (name, sequence)
#     """

#     with open(sequence) as f:
#         lines = f.readlines()

#     sequences = [(lines[i][1:].strip(), lines[i+1].strip()) \
#                     for i in range(len(lines)) if lines[i].startswith('>')]

#     return sequences


def is_protein(sequence):
    """
    Returns True if the provided sequence is a protein
    """
    amino_acids = list('MPWRLHIDEVNQSKFY')

    for aa in amino_acids:
        if aa in sequence:
            return True

    return False


def get_originals(sequences):
    """
    Make DnaOptimizationProblems for all the original sequences
    
    Input
    -----
    sequences : list
        List with SeqRecord objects with the original sequences, can be
        either DNA or protein

    Output
    -----
    originals : list
        List with SeqRecord objects with the DnaOptimizationProblems evaluations
        (score, constraints_text, objectives_text) added as annotations if they
        are DNA sequences, else, the list contains None elements
    
    """

    originals = []

    for sequence in sequences:

        name = sequence.name
        seq = str(sequence.seq)

        if not is_protein(seq):
            print("DNA sequence identified...\n")

            problem = run(seq, opt=False)

            sequence.annotations["Objectives score"] = \
                problem.objectives_evaluations().scores_sum()

            sequence.annotations["Constraints text summary"] = \
                problem.constraints_text_summary()

            sequence.annotations["Objectives text summary"] = \
                problem.objectives_text_summary()

            originals.append(sequence)
        
        else:
            print("Protein sequence identified...\n")
            originals.append(None)

    return originals


def run(sequence, opt=True):
    """
    Calls the relevant dnachisel methods to perform codon optimization

    Input
    -----
    sequence : str
        DNA sequence to be optimized
    opt : boolean
        Wether to optimize the sequence or only generate the problem. This is
        in case you only want the evaluations of the sequence
        (e.g. for the original)

    Output
    -----
    problem.sequence : str
        Optimized sequence
    """

    # Define the optimization problem

    problem = dn.DnaOptimizationProblem(
        sequence = sequence,
        constraints = [spec.AvoidPattern("GGRGG"), spec.AvoidPattern("GGTCTC"),
                    spec.AvoidPattern("GAGACC"), spec.AvoidPattern("GCGATG"),
                    spec.AvoidPattern("CATCGC"), spec.AvoidPattern("GCTCTTC"),
                    spec.AvoidPattern("GAAGAGC"), spec.AvoidPattern("CGTCTC"),
                    spec.AvoidPattern("GAGACG"), spec.AvoidPattern("GAAGAC"),
                    spec.AvoidPattern("GTCTTC"),
                    spec.EnforceGCContent(0.46,0.54),
                    spec.UniquifyAllKmers(k=8),
                    spec.EnforceTranslation()],
        objectives = [spec.CodonOptimize(species='e_coli', 
            method='match_codon_usage')]
        )

    if opt:
        # Solve the constraints, optimize with respect to the objective
        # problem.max_random_iters = 20000
        try:
            problem.resolve_constraints()
            problem.optimize()
        except dn.NoSolutionError as error:
            print('>>>>>>   NoSolutionError in one of the sequences')
            print('')
            print(error)
            print('')

    return problem


# def save_all(og_sequence, opt_sequences, destination, name):
#     """
#     Function to save the results for all the produced sequences for comparison

#     Input
#     -----
#     og_sequence : tuple
#         Original sequence, in a tuple with
#         (sequence, score, constraints_text, objectives_text)
#     opt_sequences : list
#         List of tuples with the optimized sequences in the form
#         (sequence, score, constraints_text, objectives_text)
#     destination : str
#         Directory path for saving the results
#     name : str
#         Name of the sequence
#     """

#     with open(os.path.join(
#             destination, name + '_summarized_optimization.txt'), 'w') as f:
        
#         f.write('# Best optimized sequence for each back translation\n\n')

#         if og_sequence:
#             f.write('> Original_sequence_score_%.4f\n' % (og_sequence[1]))
#             f.write(og_sequence[0] + '\n')
#             f.write('\n')
        
#         for i in range(len(opt_sequences)):
#             f.write('> %03d_score_%.4f\n' % (i+1, opt_sequences[i][1]))
#             f.write(opt_sequences[i][0] + '\n')
#             f.write('\n')

#     with open(os.path.join(
#             destination, name + '_full_optimization.txt'), 'w') as f:
        
#         f.write('# Best optimized sequence for each back translation\n\n')

#         if og_sequence:
#             f.write('> Original_sequence\n')
#             f.write(og_sequence[0] + '\n')
#             f.write(og_sequence[2] + '\n')
#             f.write(og_sequence[3] + '\n')

#         for i in range(len(opt_sequences)):
#             f.write('> %03d\n' % (i+1))
#             f.write(opt_sequences[i][0] + '\n')
#             f.write(opt_sequences[i][2] + '\n')
#             f.write(opt_sequences[i][3] + '\n')


def save_fasta(optimized, destination, jobname):
    """
    Save the optimized sequence to a file

    Input
    -----
    optimized : list
        List with tuples (name, best_sequence), where best_sequence is a tuple
        of the form (sequence, score, constraints_text, objectives_text)
    destination : str
        Directory path for saving the results
    jobname : str
        Name of the job for the results file
    """

    with open(os.path.join(destination, jobname+'_optimized.txt'), 'w') as f:
        for name, data in optimized:
            f.write('>%s_score_%.4f\n' % (name, data[1]))
            f.write(data[0]+'\n')
            f.write('\n')


def detect_ATA_ATA(sequence):
    """
    Detects the presence of two consecutive ATA codons in the given DNA sequence

    Input
    -----
    sequence : str
    """

    codon_positions = list(range(0,len(sequence), 3))
    consecutive_ATA = 0
    
    for index in codon_positions:
        icodon = sequence[index:index+3]
        
        if icodon == 'ATA':

            if index+3 < len(sequence):
                jcodon = sequence[index+3:index+6]
            
                if jcodon == 'ATA':
                    return True

    return False


def call_dnachisel(seq, n):
    """
    Calls dnachisel to optimize the provided sequence

    Input
    -----
    seq : str
        DNA sequence to be optimized
    n : int
        Number of times to run the optimization

    Output
    -----
    The information from the best optimized sequence as a tuple
    (sequence, score, constraints_text, objectives_text)
    """

    opt_sequences = []

    for j in range(n):
        problem = run(seq)
        
        score = problem.objectives_evaluations().scores_sum()

        opt_sequences.append((problem.sequence, score,\
            problem.constraints_text_summary(),\
            problem.objectives_text_summary()))

    opt_sequences = sorted(opt_sequences, key=itemgetter(1), 
        reverse=True)

    for opt in opt_sequences:
        if not detect_ATA_ATA(opt[0]):
            return opt

    return call_dnachisel(seq, n)


def optimize(in_sequence, intype, nback, n, seqid):
    """
    Executes the optimization problem n number of times, and saves the results
    
    Input
    -----
    in_sequence : str
        dna sequence to be optimized, or name of text or fasta file with
        fasta sequences
    intype : str
        If `sequence` is a file, indicate the file format as 'fasta' or 
        'genbank'
    nback : int
        number of random back-translated sequences to produce from the original
    n : int
        number of optimized sequences to produce
    seqid : str
        ID for the sequence, if it is a sequence provided as a string argument
        when calling the script, i.e. not from a file.

    Output
    -----
    originals : list
        List with SeqRecord objects with the original sequences, plus
        annotations for the dnachisel score and summaries if they are DNA
        sequences. If the original sequences are proteins, this list contains
        None objects
    optimized : list
        List with deep copies of the SeqRecord objects from the original
        sequences. The .seq attribute is changed to the optimized sequence, and
        annotations for the dnachisel score and summaries are added
    """
    
    if intype:
        sequences = list(SeqIO.parse(in_sequence, intype))
    elif os.path.isfile(in_sequence):
        raise TypeError('''Please indicate the file type with the --intype 
            argument''')
    else:
        sequences = [SeqRecord(Seq(in_sequence), id=seqid, description='')]

    # List with SeqRecords for the original sequences and their evaluations in
    # the 'annotations' dictionary
    # if they are DNA sequences, else, None
    originals = get_originals(sequences)

    # List with (SeqRecord, [backtranslated1, backtranslated2, ...])
    # for each input sequence
    backtranslated = back_translate(sequences, nback)

    # List with SeqRecords for each input sequence
    optimized = []

    for i in range(len(backtranslated)):

        sequence = backtranslated[i][0]
        btr_seqs = backtranslated[i][1]

        # List with tuples 
        # (optimized_sequence, score, constraints_text, objectives_text)
        optimized_backtranslated = [call_dnachisel(seq, n) for seq in btr_seqs]

        optimized_backtranslated = sorted(optimized_backtranslated,
            key=itemgetter(1), reverse=True)

        # Deprecated
        # if saveall:
        #     save_all(originals[i], optimized_backtranslated, destination, 
        #         sequence.name)

        best_optimized = optimized_backtranslated[0]

        sequence.seq = Seq(best_optimized[0], 
            IUPAC.unambiguous_dna)

        sequence.annotations["Objectives score"] = best_optimized[1]
        sequence.annotations["Constraints text summary"] = best_optimized[2]
        sequence.annotations["Objectives text summary"] = best_optimized[3]
        
        optimized.append(sequence)

    # Turn these methods into tests later
    assert verify_translation(originals, optimized)
    assert verify_notequal(originals, optimized)

    return originals, optimized


# def get_standard_table():
#     """
#     Gets the standard codon table from Bio.Data.CodonTable

#     Output
#     -----
#     standard_codons : dict
#         Dictionary of the form {'F':['TTT', 'TTC', ...], ...}
#     """

#     # Dictionary of the form {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', ...}
#     standard_table=CodonTable.unambiguous_dna_by_name["Standard"].forward_table

#     # Invert the standard table dictionary to make the keys the amino acids, and
#     # the values lists of possible codons
#     standard_codons = {}
#     for key, value in standard_table.items():
#         if value not in standard_codons:
#             standard_codons[value] = [key]
#         else:
#             standard_codons[value].append(key)

#     return standard_codons


def get_codon_probs():
    """
    Get dictionary with codon probabilities, in the form:

    {'*': [('TAA', 'TAG', 'TGA'), (0.64, 0.07, 0.29)], 'A': ... }
    
    """
    # Dictionary of the form
    # {'*': {'TAA': 0.64, 'TAG': 0.07, 'TGA': 0.29}, 'A': ... }
    codon_table = pct.get_codons_table("e_coli")

    # A modification is necessary since the probabilities of G do not
    # sum 1 ... substracted 0.01 to the value
    codon_table['G']['GGT'] = 0.33

    # Change it to
    # {'*': [('TAA', 'TAG', 'TGA'), (0.64, 0.07, 0.29)], 'A': ... }
    codons_probabilities = {}
    for aa, codons_ps in codon_table.items():
        codons, ps = zip(*codons_ps.items())
        codons_probabilities[aa] = [codons, ps]

    return codons_probabilities


# def random_backtranslate(seq, codon_table):
#     """
#     Performs a back translation with random codons

#     Input
#     -----
#     seq : str
#         The protein sequence to be back translated

#     codon_table : dict
#         Dicitionary in the form {'F':['TTT', 'TTC', ...], ...}
#     """

#     bt_seq = ''
#     for aa in seq:
#         codons = codon_table[aa]
#         bt_seq += np.random.choice(codons)

#     return bt_seq


def opt_backtranslate(seq, codons_probabilities, mode):
    """
    Performs a back translation taking into account the codon frequency of a
    reference organism

    Input
    -----
    seq : str
        The protein sequence to be back translated

    codons_probabilities : dict
        The probabilities for each codon for each amino acid, in the form
        {'*': [('TAA', 'TAG', 'TGA'), (0.64, 0.07, 0.29)], 'A': ... }

    mode : string
        The mode in which the sequence will be back-translated ('optimize' or 
        'random')

    Output
    -----
    bt_seq : string
        back-translated amino acid sequence
    """

    bt_seq = ''

    if mode == 'optimize':
        for aa in seq:
            codons = codons_probabilities[aa][0]
            ps = codons_probabilities[aa][1]
            bt_seq += np.random.choice(codons, p=ps)
    
    elif mode == 'random':
        for aa in seq:
            codons = codons_probabilities[aa][0]
            bt_seq += np.random.choice(codons)

    return bt_seq


def back_translate(sequences, nback, mode='optimize'):
    """
    Provides a list of randomly back-translated sequences 
    
    Input
    -----
    sequences : list
        List with SeqRecord objects for the original sequences
    nback : int
        Number of back-translated sequences to produce
    mode : str
        Mode of back translation. If 'optimized', the back translation will
        assign codons according to the codon frequency table for E. coli. 
        If 'random', the codons will be assigned randomly from the standard
        table. 

    Output:
    -----
    backtranslated : list
        List with 
        (SeqRecord, [backtranslated1, backtranslated2, backtranslated3, ...])
        tuples for each sequence
    """

    backtranslated = []

    codons_probabilities = get_codon_probs()

    for og_sequence in sequences:

        sequence = copy.deepcopy(og_sequence)

        seq = str(sequence.seq)

        if not is_protein(seq):
            bt_seqs = [seq]             # The first is the original
            seq = translate(seq)
        else:
            bt_seqs = []

        bt_seqs += [opt_backtranslate(seq, codons_probabilities, mode) \
                        for i in range(nback-len(bt_seqs))]

        backtranslated.append((sequence, bt_seqs))

    return backtranslated


def verify_translation(originals, optimized):
    """
    Verify that the optimized sequences code for the same protein as the
    originals
    """
    if all(originals):
        for i in range(len(originals)):
            if not translate(originals[i].seq) == translate(optimized[i].seq):
                return False
        return True
    else:
        return True

def verify_notequal(originals, optimized):
    """
    Verify that the optimized sequences are not the same as the originals
    """
    if all(originals):
        for i in range(len(originals)):
            if originals[i].seq == optimized[i].seq:
                return False
        return True
    else:
        return True