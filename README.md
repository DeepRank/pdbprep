Prepare PDB files for MD minimization with OpenMM Amber14 forcefield.

# 1. Install dependencies

Install the following dependencies in the same Python environment:

1. Install OpenMM: `conda install -c conda-forge openmm`
1. Install chardet: `conda install chardet`
1. Install pdb2pqr: `pip install pdb2pqr`
1. Install @joaomcteixeira fork of the `Pras_Server` as follows:

```bash
# ensure to install the fork outside the `pdbprep` repository, if needed,
# navigate to the parent directy in case you are in the `pdbprep` repository
# folder
cd ..

# clone the fork and compile the software
git clone https://github.com/joaomcteixeira/Pras_Server
cd Pras_Server
git fetch
git checkout nolog
cd Pras_Server_C++
g++ -std=c++17 src/*.cpp -o PRAS
```

Why the fork? Because in the fork's branch the logging operations were removed
to avoid writing thousand of log files to disk. All credit about PRAS should be
given to the original authors:

* https://github.com/osita-sunday-nnyigide/Pras_Server
* https://pubs.acs.org/doi/10.1021/acs.jcim.2c00571

# 2. Give execution permission to files

In the current `pdbprep` repository folder, give the necessary permissions to
the `*.py` and `.sh` files:

- If you're still on the Pras folder: `cd ../pdbprep`
- Give permission: `chmod u+x *.py pdb_prepare.sh`

The `pdb_*.py` files were adapted from the `pdb-tools` project from Alexandre
Bonvin lab.

http://www.bonvinlab.org/pdb-tools/

Here, @joaomcteixeira modifed the scripts reducing their versatility to
improve their speed. Hence, the `pdb-tools` scripts provided here won't work
outside the `pdbprep` context. If you want to use `pdb-tools` for any other need,
install the official package `pip install pdb-tools`, and cite the original work.

Add the absolute path to the PRAS file in the `pdb_prepare.sh` file (edit line 4 of the
`.sh` file).

# 3. Prepare PDBs

## 3.1 source the `setup.sh` file

From within the `pdbprep` folder, source the `setup.sh` file: `source setup.sh`.
You need to perform this operation every time you want to use `pdbprep` in a new
terminal window.

## 3.2 Prepare the PDB files

To prepare the PDB files:

1. Navigate to the folder where you want the new PDBs to be saved.
1. create a file with the list of paths to the input PDB files.
You can use `ls path/to/my/pdbs/*.pdb > pdblist` from command-line
or create a new file and paste `path/to/my/pdbs/*.pdb` inside.
1. run `pdb_prepare.sh pdblist <N>`, where `N` is the numbers of threads you want to use.

The script will create a series of folders for the different steps. If
everything goes okay, temporary PDBs will be deleted and only those in the last
folder `4_ready_to_minimize` will be saved. If something goes wrong with a
PDB, its intermediate temporary files won't be deleted and we can check them
afterwards.

Several comments in the `pdb_prepare.sh` file explain the process.


# 4. Troubleshooting

## 4.1 Can't find Python

In case the script can find the Python interpreter, type `whereis python` (or
`which python`) and update the python path in the `shabang` of the `*.py` files
accordingly.