"""
Optimizer class that will contain the methods for the codon optimization    
"""
import copy
from operator import itemgetter
from sequence_utils import back_translate, check_nucleotide, detect_ATA_ATA
from logger_utils import get_logger

import dnachisel
import dnachisel.builtin_specifications as spec
from dnachisel.DnaOptimizationProblem.DnaOptimizationProblem import DnaOptimizationProblem

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, translate


logger = get_logger(__name__)


class Optimizer:
    """
    Class that will contain the methods, constraints and goals for each
    optimization problem
    """
    
    CONSTRAINTS = [
        spec.AvoidPattern("GGRGG"), spec.AvoidPattern("GGTCTC"),
        spec.AvoidPattern("GAGACC"), spec.AvoidPattern("GCGATG"),
        spec.AvoidPattern("CATCGC"), spec.AvoidPattern("GCTCTTC"),
        spec.AvoidPattern("GAAGAGC"), spec.AvoidPattern("CGTCTC"),
        spec.AvoidPattern("GAGACG"), spec.AvoidPattern("GAAGAC"),
        spec.AvoidPattern("GTCTTC"),
        
        # promoter patterns: http://parts.igem.org/Help:Promoters/Prokaryotic_RNAP
        spec.AvoidPattern('TTGACA.{17}TATAAT'),
        spec.AvoidPattern('NNNNMRNRYTGGCACGNNNNTTGCWNNWNNNNN'),
        spec.AvoidPattern('NTCNCCCTTGAA.{17}CCCCATTTA'),
        spec.AvoidPattern('TAAAGWWY{11,12}RYCGAWRN'),
        spec.AvoidPattern('TGGATAAACATTTCACCACTGTAAGGAAAATAATTCTTATTTCGATTGTCCTTTTTACCCT'),
        
        spec.AvoidPattern('TAAGAG'), # strong RBS (Twist)
        spec.AvoidPattern('TTTTT'), spec.AvoidPattern('AAAAA'), # terminator (Twist)
        spec.AvoidPattern('10xA'), spec.AvoidPattern('10xT'),   # homopolymers >9 (Twist)
        spec.AvoidPattern('10xG'), spec.AvoidPattern('10xC'),
        spec.UniquifyAllKmers(k=20), # repeats > 20 (Twist)
        spec.EnforceGCContent(0.38,0.62), # Twist recommendation 35 - 65%

        spec.EnforceTranslation()
    ]

    OBJECTIVES = [
        spec.CodonOptimize(species='e_coli', method='match_codon_usage'),
        spec.UniquifyAllKmers(k=8, boost=5.),
        spec.EnforceGCContent(0.45,0.55, boost=3.)
    ]

    MAX_RANDOM_ITERS = 5000
    BACK_TRANSLATION_MODE = 'optimize'
    NUMBER_BACK_TRANSLATED = 2
    NUMBER_OPTIMIZATIONS = 2


    def __init__(self, sequence: SeqRecord):
        """
        Generate a DnaOptimizationProblem for a given sequence

        Args:
            sequence (SeqRecord): SeqRecord of the original sequence to be
                optimized
        """
        self.input_SeqRec = sequence
        self.sequence = str(self.input_SeqRec.seq)
        self.is_nucleotide = check_nucleotide(self.sequence)

        if self.is_nucleotide:
            logger.info('DNA sequence identified...')
            
            self.protein_sequence = translate(self.sequence)
            self.backtranslated_sequences = [self.sequence]
            self._initial_problem = self._define_DnaOptProblem(self.sequence)
            self._annotate_sequence(self.input_SeqRec, self._initial_problem)

            logger.info(f'Initial score: {self.input_SeqRec.annotations["Objectives score"]}')

        else:
            logger.info('Protein sequence identified...')
            self.protein_sequence = self.sequence
            self.backtranslated_sequences = []

        logger.info(f'Creating {self.NUMBER_BACK_TRANSLATED} back-translations...')

        self.backtranslated_sequences += [
            back_translate(self.protein_sequence, self.BACK_TRANSLATION_MODE) for i in \
                range(self.NUMBER_BACK_TRANSLATED - len(self.backtranslated_sequences))]


    def _define_DnaOptProblem(self, sequence:str) -> DnaOptimizationProblem:
        """
        Defines the DnaOptimizationProblem for a given sequence

        Args:
            sequence (str): DNA sequence to be optimized

        Returns:
            DnaOptimizationProblem: Initialized problem
        """
        problem = dnachisel.DnaOptimizationProblem(
            sequence = sequence,
            constraints = self.CONSTRAINTS,
            objectives = self.OBJECTIVES)
        
        problem.max_random_iters = self.MAX_RANDOM_ITERS

        return problem

    def _optimize_DnaOptProblem(self, problem:DnaOptimizationProblem
        ) -> DnaOptimizationProblem:
        """
        Optimize the DnaOptimizationProblem according to the goals and constraints

        Args:
            problem (DnaOptimizationProblem): Problem to optimize

        Returns:
            DnaOptimizationProblem: Optimized problem
        """
        try:
            problem.resolve_constraints()
            problem.optimize()
        except dnachisel.NoSolutionError as error:
            logger.error('NoSolutionError in one of the sequences')
            logger.error(error)

        return problem

    @staticmethod
    def _annotate_sequence(sequence:SeqRecord,
        problem:DnaOptimizationProblem):
        """
        Annotates a SeqRecord with the scores from the optimization problem

        Args:
            sequence (SeqRecord): SeqRecord to be annotated
            problem (DnaOptimizationProblem): Problem that contains the scores
                and text summaries
        """
        sequence.annotations["Objectives score"] = (
            problem.objectives_evaluations().scores_sum())
        sequence.annotations["Constraints text summary"] = (
            problem.constraints_text_summary())
        sequence.annotations["Objectives text summary"] = (
            problem.objectives_text_summary())

    def _optimize_single_sequence(self, sequence:str) -> tuple:
        """
        Optimize the DnaOptProblem for a single sequence, `NUMBER_OPTIMIZATIONS`
        times

        Args:
            sequence (str): Sequence to be optimized

        Returns:
            tuple: tuple of the form (score, problem) for the optimized sequence
                with the best score.
        """
        problem = self._define_DnaOptProblem(sequence)
        optimized_problems = []

        for i in range(self.NUMBER_OPTIMIZATIONS):
            iproblem = copy.deepcopy(problem)
            iproblem = self._optimize_DnaOptProblem(iproblem)
            optimized_problems.append(
                (iproblem.objectives_evaluations().scores_sum(), iproblem))
        
        optimized_problems = sorted(optimized_problems, key=itemgetter(0),
            reverse=True)
        
        for opt in optimized_problems:
            if not detect_ATA_ATA(opt[1].sequence):
                return opt
        
        return self._optimize_single_sequence(sequence)

    def optimize(self) -> SeqRecord:
        """
        Performs optimizations for all the back translated sequences, and
        returns the sequence with the best score

        Returns:
            SeqRecord: Optimized sequence with the best score from all the back
                translations
        """
        logger.info(f'Performing {self.NUMBER_OPTIMIZATIONS} optimizations per back-translated sequence...')
        
        optimized_backtranslated = [self._optimize_single_sequence(seq) \
            for seq in self.backtranslated_sequences]

        optimized_backtranslated = sorted(optimized_backtranslated,
            key=itemgetter(0), reverse=True)
        
        best_optimized_problem = optimized_backtranslated[0][1]

        self.best_optimized_SeqRec = copy.deepcopy(self.input_SeqRec)
        self.best_optimized_SeqRec.seq = Seq(best_optimized_problem.sequence)
        self.best_optimized_SeqRec.annotations["molecule_type"] = "DNA"
        self._annotate_sequence(self.best_optimized_SeqRec, best_optimized_problem)

        logger.info(f'Optimized sequence score: {self.best_optimized_SeqRec.annotations["Objectives score"]}')

        return self.best_optimized_SeqRec

