#!/usr/bin/env python
"""
Add hydrogens to models
"""
import os
import sys
import json
from copy import copy
from pathlib import Path

import openmm as mm
import openmm.app as app
import openmm.unit as units


def write_structure(fout_name, structure):
    with open(fout_name, 'w') as fout:
        app.PDBFile.writeFile(structure.topology, structure.positions, fout)


# PARAMETERS
forcefield_model = 'amber14-all.xml' #'charmm36.xml'
water_model = 'amber14/tip3p.xml' #'charmm36/tip3p-pme-b.xml'
kcal_mole = units.kilocalorie_per_mole

platform = mm.Platform.getPlatform(1)
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

structure.positions = model.positions
structure.topology = model.topology

#
# Setup simulation
#
system = forcefield.createSystem(structure.topology)

integrator = mm.LangevinIntegrator(
    310*units.kelvin,
    1.0/units.picosecond,
    2.0*units.femtosecond,
    )

integrator.setRandomNumberSeed(917)
integrator.setConstraintTolerance(0.00001)

simulation = app.Simulation(structure.topology, system, integrator)

context = simulation.context
context.setPositions(model.positions)

state = context.getState(getEnergy=True)
ini_ene = state.getPotentialEnergy().value_in_unit(kcal_mole)

simulation.minimizeEnergy(maxIterations=100)
structure.positions = context.getState(getPositions=True).getPositions()

state = context.getState(getEnergy=True)
#min_ene = state.getPotentialEnergy().value_in_unit(kcal_mole)

#
# Minimize
#

state = context.getState(getEnergy=True)
ini_ene = state.getPotentialEnergy().value_in_unit(kcal_mole)

simulation.minimizeEnergy(maxIterations=100)
structure.positions = context.getState(getPositions=True).getPositions()

state = context.getState(getEnergy=True)
#min_ene = state.getPotentialEnergy().value_in_unit(kcal_mole)


fout_fname = Path(output_folder, Path(input_structure).name)
write_structure(fout_fname , structure)
