# CoolerCodonOpt
---
Codon optimization using the [dnachisel](https://pypi.org/project/dnachisel/) library

## Installation with conda (recommended)

Conda is an open source package and environment management system that runs on Windows, macOS and Linux. To know more about Conda take a look at the workshop [Introduction to Conda for (Data) Scientists](https://carpentries-incubator.github.io/introduction-to-conda-for-data-scientists/) or the [YouTube workshop](https://www.youtube.com/channel/UCR1RFwgvADo5CutK0LnZRrw/featured) prepared by the KAUST Visualization Core Lab.

1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) in your system:

   For MacOS:

   -  Download the installer in your home directory: https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

   - Open your terminal, go to your home directory (or wherever you downloaded the file) and run the installation script:

     ```
     bash Miniconda3-latest-MacOSX-x86_64.sh
     ```

   - The script will present several prompts that allow you to customize the Miniconda install. It is recommended to accept the default settings. However, when prompted with the following:

     ```
     Do you wish the installer to initialize Miniconda3
     by running conda init?
     ```

     Type `yes` to avoid to manually initialize Conda later. If you accidentally type `no`, when the script finishes you simply need to type the following command:

     ```
     conda init bash
     ```

   - Close and open the terminal again. You will see either `(miniconda3)` or `(base)` before the terminal prompt, which means conda was installed correctly and the `base` environment is activated.


2. Open your terminal and go to the directory where you would like to save the project. Then run the following command: 

   `git clone https://github.com/strubelab/CoolerCodonOpt.git`
   
   - The directory `CoolerCodonOpt/` is created.

3. Create a virtual environment to work on:

   - Go to the `CoolerCodonOpt` directory
   
     `cd CoolerCodonOpt`
   
   - Create a virtual environment
      
     `conda env create --prefix ./env --file environment.yml`
   
   - Activate the virtual environment
   
     `conda activate ./env`

4. Install the `CoolerCodonOpt` package in your virtual environment:

    ```
    pip install -e .
    ```     

## Installation with venv and pip

1. Verify that you have python 3 installed in your computer:
   
   - Open your terminal and type `python --version`
   - If the version number is less than 3, for example `Python 2.7.10`, download the current python version from https://www.python.org/ or alternatively follow the installation instructions in https://realpython.com/installing-python/
   
   
2. Open your terminal and go to the directory where you would like to save the project. Then run the following command: 

   `git clone https://github.com/strubelab/CoolerCodonOpt.git`
   
   - The directory `CoolerCodonOpt/` is created.
   
   
3. Create a virtual environment to work on:

   - Go to the `CoolerCodonOpt` directory
   
     `cd CoolerCodonOpt`
   
   - Create a virtual environment
      
     `python3 -m venv venv`
   
   - Activate the virtual environment
   
     `source venv/bin/activate`
     
4. Install requirements and the `CoolerCodonOpt` package

    ```
    pip install -r requirements.txt
    pip install -e .
    ```
   
5. When you finish working, you can deactivate the virtual environment
      
     `deactivate`

## Usage

1. Now you can call the `optimize` script from any location:

```
optimize --help
```
   
```
usage: optimize [-h] [-v] [-d DESTINATION] input

Takes a DNA sequence and optimizes it for expression in E. coli

positional arguments:
  input                 Fasta file with the sequence(s) to be optimized, or a directory with fasta files.

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Show the constraints evaluations, and optimization objectives score.
  -d DESTINATION, --destination DESTINATION
                        Path for saving the resulting sequences. It defaults to the same directory as the input.
```

### A. Run providing a fasta file with one or more sequences

The following command will optimize the given sequence, output the optimized sequence and optimization score to the terminal and save the optimized sequence to a file in fasta format:

``` bash
python optimize.py sequence.fasta --verbose
```

### B. Provide a directory name with fasta sequence files

The following command will take all the sequences in the `sequences/` directory, and save them in the `optimized_sequences/` directory.

``` bash
python optimize.py sequences/ --destination optimized_sequences/
```
