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
usage: optimize.py [-h] [-i INTYPE] [-o OUTTYPE] [-v] [-p] [-d DESTINATION]
                   [-q SEQID] [-j JOBNAME] [-n NUMBER] [-b NBACK]
                   sequence

Takes a DNA sequence and optimizes it for expression in E. coli

positional arguments:
  sequence              DNA sequence to be optimized. Takes a single sequence
                        as a string, or the name of a file with one or more
                        sequences in fasta or genbank format. If a file is
                        provided, the --intype argument has to be given as
                        well.

optional arguments:
  -h, --help            show this help message and exit
  -i INTYPE, --intype INTYPE
                        If `sequence` is a file, indicate the file type as
                        `fasta` or `genbank`.
  -o OUTTYPE, --outtype OUTTYPE
                        Save the resulting optimized sequences to a file, in
                        the indicated format. It can be `fasta`, `genbank` or
                        `fasta/genbank` for both.
  -v, --verbose         Show the constraints evaluations, and optimization
                        objectives score.
  -p, --separate        Save the optimized sequences from the input file in
                        separate files, as opposed to all in one single file
  -d DESTINATION, --destination DESTINATION
                        Path for saving the resulting sequences. --outtype has
                        to be provided.
  -q SEQID, --seqid SEQID
                        Name/identifier of the sequence. Only used when a
                        single sequence is provided as a string, and you want
                        to save it to a file.
  -j JOBNAME, --jobname JOBNAME
                        Name of the job that will go in the name of the
                        results file. --outtype has to be provided.
  -n NUMBER, --number NUMBER
                        Number of optimized sequences to produce per each
                        back-translated one (default=10).
  -b NBACK, --nback NBACK
                        Number of random back-transalted sequences to be
                        produced (default=10). 0 to not back translate and
                        optimize only from the original sequence.
```

#### A. Run providing a single sequence as a string

The following command will optimize the given sequence, output the optimized sequence and optimization score to the terminal and save the optimized sequence to file in fasta format:

``` bash
python optimize.py ASFGTFTGFAWDWDWEERFFKLKLIP -v --outtype fasta
```

#### B. Provide one or more sequences in a FASTA file

The following command will take the sequences given in `dna_sequences.fasta`, optimize them and save them to a file in genbank format in the folder `../example_data`:

``` bash
python optimize.py ../example_data/dna_sequences.fasta --intype fasta --outtype genbank --destination ../example_data
```
