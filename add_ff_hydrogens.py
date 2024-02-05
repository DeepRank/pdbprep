#!/usr/bin/env python
"""
Add hydrogens to models and perform a short minimization in vacuum.
"""
import argparse
import json
from copy import copy
from pathlib import Path

import openmm as mm
import openmm.app as app
import openmm.unit as units


def write_structure(fout_name, structure):
    with open(fout_name, 'w') as fout:
        app.PDBFile.writeFile(structure.topology, structure.positions, fout)


ap = argparse.ArgumentParser()
ap.add_argument(
    '-s',
    '--structure',
    help="Input Structure",
    required=True,
    )

ap.add_argument(
    '-p',
    '--protonation',
    help="Protonation state file.",
    required=True,
    )

ap.add_argument(
    '-o',
    '--output',
    help="Output folder.",
    required=True,
    )

ap.add_argument(
    '-m',
    '--minimize',
    help="Whether or not to perform a short minimization in vacuum.",
    action='store_true',
    )

ap.add_argument(
    '-t',
    '--temperature',
    help="Minimization temperature in Kelvin. Default 310K.",
    default=310,
    type=int,
    )

ap.add_argument(
    '-rs',
    '--random-seed',
    help='Randon seed for minimization.',
    default=917,
    type=int,
    )

ap.add_argument(
    '-it',
    '--max-iterations',
    dest='max_iterations',
    help='Maximum iterations for minimization.',
    default=100,
    type=int,
    )

cmd = ap.parse_args()

# PARAMETERS
forcefield_model = 'amber14-all.xml' #'charmm36.xml'
water_model = 'amber14/tip3p.xml' #'charmm36/tip3p-pme-b.xml'
kcal_mole = units.kilocalorie_per_mole
platform_properties = {'Threads': str(1)}

input_structure = cmd.structure
protonated_sequence_fname = cmd.protonation
output_folder = cmd.output

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

if cmd.minimize:
    system = forcefield.createSystem(structure.topology)

    integrator = mm.LangevinIntegrator(
        cmd.temperature*units.kelvin,
        1.0/units.picosecond,
        2.0*units.femtosecond,
        )

    integrator.setRandomNumberSeed(cmd.random_seed)
    integrator.setConstraintTolerance(0.00001)

    simulation = app.Simulation(
        structure.topology,
        system,
        integrator,
        platformProperties=platform_properties,
        )

    context = simulation.context
    context.setPositions(model.positions)

    state = context.getState(getEnergy=True)
    ini_ene = state.getPotentialEnergy().value_in_unit(kcal_mole)

    simulation.minimizeEnergy(maxIterations=cmd.max_iterations)
    structure.positions = context.getState(getPositions=True).getPositions()

    state = context.getState(getEnergy=True)

    simulation.minimizeEnergy(maxIterations=100)

    structure.positions = context.getState(getPositions=True).getPositions()

fout_fname = Path(output_folder, Path(input_structure).name)
write_structure(fout_fname , structure)
