Prepare PDB files for MD minimization with OpenMM Amber14 forcefield.

# 1. Install dependencies

Install the following dependencies in the same Python environment:

1. Install OpenMM: `conda install -c conda-forge openmm`
1. Install chardet and scipy: `conda install chardet scipy` 
1. Install pdb2pqr: `pip install pdb2pqr`
1. Install @joaomcteixeira fork of the `Pras_Server` as follows:
  - make sure to do this outside of the current repo (so do `cd ..` if you are currently in the folder for this repo)

```bash
git clone https://github.com/joaomcteixeira/Pras_Server
cd Pras_Server
git fetch
git checkout nolog
cd Pras_Server_C++
g++ -std=c++17 src/*.cpp -o PRAS
```

Why the fork? Because in the fork's branch the logging was removed to avoid
writing thousand of log files to disk. All credit about PRAS should be given to
the original authors:

* https://github.com/osita-sunday-nnyigide/Pras_Server
* https://pubs.acs.org/doi/10.1021/acs.jcim.2c00571

PRAS performs better to add missing heavy atoms than other more known programs.

# 2. Give permission
In the current repo (i.e. not the Pras repo),  give the necessary permissions to the `*.py` and `.sh` files.

- If you're still on the Pras folder: `cd ../pdbprep` 
- Give permission: `chmod u+x *.py pdb_prepare.sh`

The `pdb_*.py` files are taken from the `pdb-tools` project from Alexandre Bonvin lab.

http://www.bonvinlab.org/pdb-tools/

Here, @joaomcteixeira modifed the scripts to reduce their versatility and
improve their speed. Hence, the `pdb-tools` scripts provided here won't work
outside this context. If you want to use `pdb-tools` for any other need,
install the official package `pip install pdb-tools`.

If PRAS was not placed in the default location (as suggested above), 
add the absolute path to the PRAS file in the `pdb_prepare.sh` file (edit its line 4).

## Add the folder with the scripts to `PATH`:

`echo "export PATH=\$PATH:`pwd`" >> .bashrc`

Source `.bashrc` after.
- Note that this step must always be done (in the correct path), each time you restart the terminal.

## Update `*.py` shebang

Type `whereis python` (or `which python`) and update the python path in the
`shabang` of the `*.py` files accordingly.

# 3. Prepare PDBs

To prepare PDBs:

1. Navigate to the folder where you want the new PDBs to be saved.
1. create a file with the list of paths to the input PDB files. You can use `ls path/to/my/pdbs/*.pdb > pdblist` from commandline or create a new file and paste `path/to/my/pdbs/*.pdb` inside.
1. run `pdb_prepare.sh pdblist <N>`, where `N` is the numbers of threads you want to use.

The script will create a series of folders for the different steps. If
everything goes okay, temporary PDBs will be deleted and only those in the last
folder `4_ready_to_minimize` will be saved. If something goes wrong with a
PDB, its intermediate temporary files won't be deleted and we can check them
afterwards.

Several comments in the `pdb_cleaning_v2.sh` file explain the process.
