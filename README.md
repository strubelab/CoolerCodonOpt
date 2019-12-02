# CoolerCodonOpt
---
Codon optimization using the [dnachisel](https://pypi.org/project/dnachisel/) library

## Installation

1. Verify that you have python 3 installed in your computer:
   
   - Open your terminal and type `python --version`
   - If the version number is less than 3, for example `Python 2.7.10`, download the current python version from https://www.python.org/ or alternatively follow the installation instructions in https://realpython.com/installing-python/
   
   
2. Open your terminal and go to the directory where you would like to save the project. Then run the following command: 

   `git clone https://github.com/strubelab/CoolerCodonOpt.git`
   
   - The directory `CoolerCodonOpt/` is created.
   
   
3. Create a virtual environment to work on:

   - Go to the `CoolerCodonOpt` directory
   
   `cd CoolerCodonOpt`
   
   - Install `virtualenv` if not present, and create the virtual environment
   
   `pip install virtualenv`
   
   `virtualenv --python=python3 venv`
   
   - Activate the virtual environment
   
   `source venv/bin/activate`
   
   - Deactivate the virtual environment (when you finish working)
   
   `deactivate`
   
   
4. Install dnachisel with the following command:

   `pip install dnachisel[reports]`
   

## Usage

1. Go to the CoolerCodonOpt/scripts directory and run the script:

   `cd CoolerCodonOpt/scripts`
   
   `python optimize.py -h`
   
   
```
usage: optimize.py [-h] [-v] [-b NBACK] [-n NUMBER] [-s] [-e] [-d DESTINATION]
                   [-a] [-q SEQNAME] [-m MODEBT] [-j JOBNAME]
                   sequence

Takes a DNA sequence and optimizes it for expression in E. coli

positional arguments:
  sequence              DNA sequence to be optimized. Takes a single sequence
                        as a string, or a text file with FASTA sequences.

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Show the constraints evaluations, and optimization
                        objectives score.
  -b NBACK, --nback NBACK
                        Number of random back-transalted sequences to be
                        produced (default=10). 0 to not back translate and
                        optimize only from the original sequence.
  -n NUMBER, --number NUMBER
                        Number of optimized sequences to produce per each
                        back-translated one (default=10).
  -s, --save            Save the resulting optimized sequences to a file (only
                        the highest scoring).
  -e, --saveall         Save the best optimized sequences from each back-
                        translation for comparison.
  -d DESTINATION, --destination DESTINATION
                        Directory path for saving the resulting sequences and
                        scores, if --save or --saveall options are selected.
  -a, --atg             Add ATG at the beginning of the sequence.
  -q SEQNAME, --seqname SEQNAME
                        Name/identifier of the sequence. Only used when a
                        single sequence is provided.
  -m MODEBT, --modebt MODEBT
                        Mode to run back translation of the sequence, it can
                        be "optimize" to take into account codon frequencies,
                        or "random" to give random codons (default=optimize).
  -j JOBNAME, --jobname JOBNAME
                        Name of the job that will go in the name of the
                        results file when activating --save.
```

#### A. Run providing a single sequence as a string

``` bash
python optimize.py ATGGGCGCTGGGCGTTTCGCCCCCCGTGTCGCCTGCACAAGGGCGACATGCTACGTGAAA -v
```

#### B. Provide one or more sequences in a FASTA file

``` bash
python optimize.py ../example_data/sequences.fasta --save --jobname example1
```
