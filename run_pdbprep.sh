#!/bin/bash

# replace this by the path of your PRAS executable
PRAS_exec="/home/joao/github/Pras_Server/Pras_Server_C++/PRAS"
pdbtools_folder="0_pdbtools"
pras_folder="1_pras"
pdb2pqr_folder="2_pdb2pqr"
protonation_folder="3_protonation"
with_Hs_folder="4_ready_to_minimize"


assert_file_exists () {
    if [ ! -x "$1" ]; then
        echo "ERROR: $1 does not exist or is not executable."
        exit 1
    fi
}


do_pdbtools() {
    # basic initial preparion of the PDBs:
    filename="${1##*/}"
    basename="${filename%.*}"
    pdbtoolsout="${pdbtools_folder}/${basename}.pdb"
    # removes non coordinate lines for simplicity
    pdb_keepcoord.py $1 |\
    # removes waters
    pdb_delresname.py |\
    # converts residue names to standard names, ex: MSE to MET
    pdb_rplresname.py |\
    # selects most probable alternative location
    pdb_selaltloc.py |\
    # fixes inserts
    pdb_fixinsert.py |\
    # sort chains and resides, necessary for OpenMM
    pdb_sort.py |\
    # renumber residues from 1
    pdb_reres.py |\
    # renumber atoms from 1
    pdb_reatom.py |\
    # tidy cleans the PDB, adds TER, etc.
    pdb_tidy.py > $pdbtoolsout
}

do_pras() {
    # Adds missing heavy atoms, which generally are many.
    # Does NOT perform any sidechain rotation or sidechain packing.
    # In this stage, metal ions and ligands are lost.
    pdbtoolsout="${pdbtools_folder}/${1}.pdb"
    prasout="${pras_folder}/${1}.pdb"
    $PRAS_exec \
        -i $pdbtoolsout \
        -o $prasout \
        -f none \
        -r yes \
        -m yes \
        -l no \
        -h no \
        -p no \
        -s no > /dev/null
}

do_pdb2pqr() {
    # Calculates the protonation states.
    pdb2pqrout="${pdb2pqr_folder}/${1}.pdb"
    pdb2pqrlog="${pdb2pqr_folder}/${1}.log"
    pdb2pqr $prasout $pdb2pqrout \
        --ff=AMBER \
        --nodebump \
        --keep-chain \
        --log-level CRITICAL > /dev/null
}

do_protonation() {
    # Reads the calculated protonation states and generates a JSON file
    # with that information. The JSON will be given to OpenMM Modeller to correctly
    # place the hydrogen atoms according to the pdb2pqr calculations.
    pdb2pqrout="${pdb2pqr_folder}/${1}.pdb"
    protonationout="${protonation_folder}/${1}.json"
    detect_protonation.py $pdb2pqrout $protonationout
}

do_add_hydrogens() {
    # Add hydrogens to the models using OpenMM Modeller and the same forcefield
    # that will be used for the MD minimization.
    prasout="${pras_folder}/${1}.pdb"
    protonationout="${protonation_folder}/${1}.json"
    add_ff_hydrogens.py \
        -s $prasout \
        -p $protonationout \
        -o $with_Hs_folder \
        -m \
        -t 310 \
        -rs 917 \
        -it 100
}

do_clean() {
    # cleans the temporary PDB files and leaves only the final PDBs with the
    # forcefield hydrogen atoms.
    pdbtoolsout="${pdbtools_folder}/${1}.pdb"
    prasout="${pras_folder}/${1}.pdb"
    pdb2pqrout="${pdb2pqr_folder}/${1}.pdb"
    pdb2pqrlog="${pdb2pqr_folder}/${1}.log"
    protonationout="${protonation_folder}/${1}.json"

    rm $pdbtoolsout $prasout $pdb2pqrout $pdb2pqrlog $protonationout
}

prepare_pdb() {

    filename="${1##*/}"
    basename="${filename%.*}"
    echo "processing ${basename}..."
    do_pdbtools $1 && \
    do_pras $basename && \
    do_pdb2pqr $basename && \
    do_protonation $basename && \
    do_add_hydrogens $basename  && \
    do_clean $basename
    # the clean step takes place only for those successfully treated PDBS.
    # PDBs where errors occurred remain in their intermediate folders.
}

if [[ -z $1 ]]; then
    echo "ERROR: You need to give an input pdb list."
    exit 1
fi
assert_file_exists $PRAS_exec

reset_folders=(\
    $pdbtools_folder \
    $pras_folder \
    $pdb2pqr_folder \
    $protonation_folder \
    )
for folder in ${reset_folders[@]}
do
    if [ -d $folder ]
    then
        rm -r $folder
    fi
    mkdir -p $folder
done

mkdir -p $with_Hs_folder

if [ -z $2 ]
then
    ncores=1
else
    ncores=$2
fi

echo "performing with ${ncores} CPUS.."
for pdb in `cat $1`
do
   if [[ -f "$with_Hs_folder/${pdb##*/}" ]]; then
       echo "${pdb##*/} exists in the $with_Hs_folder/ folder, skipping..."
   else
   prepare_pdb $pdb &
   if [[ $(jobs -r -p | wc -l) -ge $ncores ]]; then wait -n; fi
   fi
done
wait
echo "all done"
