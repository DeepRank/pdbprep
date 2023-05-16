#!/usr/bin/env python
"""
Add hydrogens to models
"""
import os
import sys
import json
from copy import copy
from pathlib import Path

import openmm
import openmm.app as app
import openmm.unit as units


def write_structure(fout_name, structure):
    with open(fout_name, 'w') as fout:
        app.PDBFile.writeFile(structure.topology, structure.positions, fout)


# PARAMETERS
forcefield_model = 'amber14-all.xml' #'charmm36.xml'
water_model = 'amber14/tip3p.xml' #'charmm36/tip3p-pme-b.xml'

platform = openmm.Platform.getPlatform(1)
threads = os.cpu_count() - 1
platform_properties = {'Threads': str(threads)}

input_structure = sys.argv[1]
protonated_sequence_fname = sys.argv[2]
output_folder = sys.argv[3]

# PREPARES MODEL
forcefield = app.ForceField(forcefield_model, water_model)
structure = app.PDBFile(input_structure)

## reads the sequence
with open(protonated_sequence_fname, 'r') as fin:
    findict = json.load(fin)
    key = list(findict.keys())[0]
    protonated_sequence = copy(findict[key])
    del findict, key

model = app.Modeller(structure.topology, structure.positions)
model.addHydrogens(forcefield=forcefield, variants=protonated_sequence)

fout_fname = Path(output_folder, Path(input_structure).name)
write_structure(fout_fname , model)
