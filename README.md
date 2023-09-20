**Prepare PDB files for MD minimization with OpenMM Amber14 forcefield.**

`pdbprep` is not a Python package. It is a series of Python scripts orchestrated
by a `bash` file that defines the preparation pipeline; hence, you can run each script
individually if needed.

# Performed tasks

`pdbpred` ensures the treated PDBs comply with the consistency standards needed
for DeepRank input. `pdbprep` is a modular pipeline that performs the following
tasks:

1. Cleans PDBs (uses scripts adapted from the
[`pdb-toos`](https://www.bonvinlab.org/pdb-tools/) package):
    1. Keeps only coordinate lines to simplify the input PDBs
    1. Removes water (`HOH`) molecules
    1. Replaces certain residue names to standard names, for example `MSE` to
    `MET` and `HIP` to `HIS`
    1. Selects the most probable alternative locations (discards others)
    1. Fixes inserts
    1. Sorts chains and residues (necessary for OpenMM)
    1. Renumber residues starting from 1
    1. Renumber atoms starting from 1
    1. Cleans the PDB (uses `pdb_tidy`)
1. Runs [PRAS](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00571) to add missing heavy atoms
1. Runs [`pdb2pqr`](https://pdb2pqr.readthedocs.io/en/latest/) to calculate the protonation states of the polar hydrogens
1. Reads the calculated protonation states and prepares a file that will serve
as input to OpenMM
1. Adds all hydrogens with OpenMM and the Amber14 forcefield using the
information provided in the previous step
1. If all steps succeeded, deletes all intermediate files and keeps only the final PDB.

# How to install

## Clone this repository

Clone this repository:

```bash
git clone https://github.com/DeepRank/pdbprep
```

## Install dependencies

`pdbpred` requires Python 3 and a `bash` shell to run. Ensure these are
installed. We suggest using [Anaconda](https://www.anaconda.com/download).

Install the following dependencies in the Python environment you wish to run
`pdbpred`.

1. Install OpenMM: `conda install -c conda-forge openmm`
1. Install chardet: `conda install chardet`
1. Install pdb2pqr: `pip install pdb2pqr`
1. Install @joaomcteixeira fork of the `Pras_Server` as follows and **outside**
the `pdbpred` folder:

```bash
# clone the fork and compile the software from the `nolog` branch
git clone https://github.com/joaomcteixeira/Pras_Server
cd Pras_Server
git fetch
git checkout nolog
cd Pras_Server_C++
g++ -std=c++17 src/*.cpp -o PRAS
```

Why installing the fork and not the original source? Because in the fork's
branch the logging operations were removed to avoid writing thousand of log
files to disk. All credit about PRAS should be given to the original authors:

* https://github.com/osita-sunday-nnyigide/Pras_Server
* https://pubs.acs.org/doi/10.1021/acs.jcim.2c00571

## Give execution permission to files

In the `pdbprep` repository folder, give the necessary permissions to the `*.py`
and `run_pdbprep.sh` files:

```bash
chmod u+x *.py run_pdbprep.sh
```

The `pdb_*.py` files were adapted from the `pdb-tools` project from Alexandre
Bonvin lab.

http://www.bonvinlab.org/pdb-tools/

Here, @joaomcteixeira modifed the scripts reducing their versatility to
improve their speed. Hence, the `pdb-tools` scripts provided here won't work
outside the `pdbprep` context. If you want to use `pdb-tools` for any other need,
install the official package `pip install pdb-tools`, and cite the original work.

If PRAS was not placed in the default location (as suggested above), add the
absolute path to the PRAS file in the `run_pdbprep.sh` file (edit line 4 of
`run_pdbprep.sh`).

# How to run

## Source the `setup.sh` file

From within the `pdbprep` folder, source the `setup.sh` file: `source setup.sh`.
You need to perform this operation every time you want to use `pdbprep` in a new
terminal window.

## Prepare a list of the target PDBs

To prepare the PDB files:

1. Navigate to the folder where you want the new PDBs to be saved.
1. create a file with the list of paths to the input PDB files. Paths can be
relative to the current folder or absolute.
You can use `ls path/to/my/pdbs/*.pdb > pdblist` to perform this operation. The
file should contain lines like the following, pointing to the input PDBs:

```
A6/6A6I.pdb
AY/7AYE.pdb
D2/7D2T.pdb
GS/6GS2.pdb
```

## Run the pipeline - prepare the PDBs!

To execute the pipeline on the list of target PDBs, run:

```
run_pdbprep.sh pdblist <N>
```

Where `N` is the number of threads (cores) you want to use. The multithreading
operation follows an *embarrassingly parallel* scheme where each thread will
take a PDB from the list and process it independently until the end.

## What to expect?

The `run_pdbprep.sh` script will create a series of numbered indexed folders to
store the temporary PDBs for the different intermediate steps (`0_*`, `1_*`,
`2_*`, ...). If the preparation succeeds, the temporary PDBs will be deleted and
only those in the last folder `4_ready_to_minimize` will be saved. If something
goes wrong with a PDB, its intermediate temporary files won't be deleted so that
errors can be traced.

At startup, `run_pdbprep.sh` will delete all temporary folders (and the files
inside), keeping only the folder with the *ready to minimize* structures.
`run_pdbpred.sh` will skip those PDBs listed in the input `pdblist` that were
already treated and are present in the *4_ready_to_minimize* folder. Therefore,
you can restart a previously halted run without needing to repeat the already
completed PDBs.

# Troubleshooting

## Can't find Python

In case the script can find the Python interpreter, type `whereis python` (or
`which python`) and update the python path in the `shabang` (1st line) of all
`*.py` files accordingly.

## Need help?

Contact us by opening a new [issue](https://github.com/DeepRank/pdbprep/issues).

# Citing

When using `pdbpred` you should acknowledge the following software, follow the
links for information about how to cite:
* [pdb-tools](https://f1000research.com/articles/7-1961)
* [Pras](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00571)
* [pbd2pqr](https://pdb2pqr.readthedocs.io/en/latest/supporting.html#citing-our-software)
* [OpenMM](https://openmm.org/)
