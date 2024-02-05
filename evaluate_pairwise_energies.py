#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate atom-pairwise vdW and electrostatic energies (kcal/mole)
between atoms in a protein complex interface using OpenMM and
CustomNonbondedForces.

If the protein complex has more than two chains, calculates the energies
for all interchain possibilities.

By default, saves to disk a `energy.log` file with the energies for all
the atom pairs plus the sum of the energies.

OpenMM-based algorithm by:
Joao Rodrigues @ Stanford, 2018
github: @JoaoRodrigues
e-mail: j.p.g.l.m.rodrigues@gmail.com

Adapted to serve a protein complex interface and CLI by:
Joao M. C. Teixeira, 2023
github: @joaomcteixeira
e-mail: joaomcteixeira@gmail.com

USAGE
-----

$ python evaluate_pairwise_energies.py -s PDB
$ python evaluate_pairwise_energies.py -h
"""
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import argparse
import itertools
import json
import logging
import operator
import os
import sys
from pathlib import Path

import openmm as mm
import openmm.app as app
import simtk.unit as units
from openmm.app.element import hydrogen as ElementHydrogen

#
# Constants
#
_zero_q = 0.0 * units.elementary_charge
_zero_sig = 0.0 * units.nanometer
_zero_eps = 0.0 * units.kilojoule_per_mole

R_IJ = 10 * units.angstrom  # cutoff for nonbonded potentials
NONB_METHOD = mm.NonbondedForce.CutoffNonPeriodic

#
# Parameter defaults
#
DEFAULT_FF = 'amber14-all.xml'
DEFAULT_ADD_HYDROGENS = False
DEFAULT_PH = 7.0
DEFAULT_MINI_SEED = 917
DEFAULT_MAX_IT = 100
DEFAULT_MINI_OUT = False
DEFAULT_CONTACT_CUTOFF = 5.0
DEFAULT_ENERGY_LOG = 'energy.log'


#
# Functions
#
def _process_hydrogen_variants(variants):
    """
    Process variants description.

    If `variants` is a JSON file, reads the residue definition list
    defined in the first key.

    If `variants` is already a list, returns that list.

    Raises `ValueError` if none of the above matches.
    """
    try:
        with open(variants, 'r') as fin:
            parsed_variants = json.load(fin)

    except (FileNotFoundError, TypeError):
        if isinstance(variants, (list, tuple)) or variants is None:
            return variants
        else:
            raise ValueError(
                'Protonation variants should be a list or a '
                'valid path to a JSON file.'
                )

    return list(parsed_variants.values())[0]


class _ProcessHydrogenVariants(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        try:
            variants_list = _process_hydrogen_variants(value)
        except ValueError as err:
            parser.error(str(err))
        setattr(namespace, self.dest, variants_list)


def write_structure(fout_name, structure):
    """Write OpenMM structure to a PDB file."""
    with open(fout_name, 'w') as fout:
        app.PDBFile.writeFile(structure.topology, structure.positions, fout)
    logging.info('Saved {}'.format(fout_name))


def sq_eucl_dij(particle_i, particle_j):
    """
    Returns the square euclidean distance between two particles
    """
    ix, iy, iz = particle_i
    jx, jy, jz = particle_j
    sq_dij = (ix-jx)**2 + (iy-jy)**2 + (iz-jz)**2
    return sq_dij.value_in_unit(units.angstrom**2)


#
# Format Logger
#
logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format='[%(asctime)s] %(message)s',
    datefmt='%Y/%m/%d %H:%M:%S',
    )

#
# Configure CLI
#
ap = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '--structure',
    '-s',
    dest='structure_fpath',
    help='Structure file',
    required=True,
    type=Path,
    )

ap.add_argument(
    '--force-field',
    '-ff',
    help='Forcefield name according to OpenMM. Defaults to Amber14-all.',
    dest='ff_name',
    default=DEFAULT_FF,
    type=str,
    )

ap.add_argument(
    '--silence-log',
    dest='silence_log',
    help=(
        'Disables logging operations to STDOUT. Useful when running '
        'highthroughput. Defaults to `False`.'
        ),
    default=False,
    action='store_true',
    )

ap_hydrogens = ap.add_argument_group(
    title='Add hydrogens',
    description=(
        'Parameters that configure adding hydrogens. '
        'If the input structure already contains hydrogens, and you want to use '
        'those, you can skip this.'
        ),
    )

ap_hydrogens.add_argument(
    '--add-hydrogens',
    dest='add_hydrogens',
    help='Whether or not to add hydrogens.',
    action='store_true',
    )

ap_hydrogens.add_argument(
    '--hydrogen-variants',
    help=(
        'The protonation states for all residues in the format of a'
        'JSON file. '
        'The JSON file contains a single key:value pair where the key can be '
        'any and the value is a list with the residue codes '
        'according to openmm.app.modeller.Modeller.addHydrogens.'
        ),
    dest='hydrogen_variants',
    type=Path,
    action=_ProcessHydrogenVariants,
    )

ap_hydrogens.add_argument(
    '--hydrogen_ph',
    dest='hydrogen_ph',
    help='The considered pH when adding hydrogens. Defaults to 7.0.',
    type=float,
    default=DEFAULT_PH,
    )

ap_minimize = ap.add_argument_group(
    title='Short vaccum minimization.',
    description='Parameters for short vacuum minimization.',
    )

ap_minimize.add_argument(
    '--minimize',
    help='Whether to perform a short minimization in vacuum.',
    action='store_true',
    )

ap_minimize.add_argument(
    '--minimization-random-seed',
    dest='minimization_random_seed',
    help='Minimization random seed.',
    default=DEFAULT_MINI_SEED,
    type=int,
    )

ap_minimize.add_argument(
    '--minimization-max-iterations',
    dest='minimization_max_iterations',
    help=(
        'The maximum numer of iterations for the short vacuum minimization. '
        'Defaults to 100.'
        ),
    default=DEFAULT_MAX_IT,
    type=int,
    )

ap_minimize.add_argument(
    '--minimized-output',
    dest='minimized_output',
    help=(
        'The path to save the minimized structure. '
        'If not given, no output will be saved.'
        ),
    default=False,
    )

ap_energy = ap.add_argument_group(
    title='Energy calculation.',
    description='Parameters that control the atom pairwise energy calcualtion.'
    )

ap_energy.add_argument(
    '--contact-cutoff',
    dest='contact_cutoff',
    help='The maximum distance to consider a contact. Defaults to 5 Angstroms.',
    default=DEFAULT_CONTACT_CUTOFF,
    type=float,
    )

ap_energy.add_argument(
    '--energy-log',
    dest='energy_log',
    help=(
        'The file to save the pairwise energies in kcal/mole. '
        'Defaults to `energy.log`. '
        'To disable writing the energy log file, provide the argument '
        'without parameters, that is, write only `--energy-log`.'
        ),
    default=DEFAULT_ENERGY_LOG,
    action='store_false',
    )


#
# Main Code
#
def calculate_interface_atom_pairwise_energies(
        structure_fpath,
        ff_name=DEFAULT_FF,
        #
        add_hydrogens=DEFAULT_ADD_HYDROGENS,
        hydrogen_variants=None,
        hydrogen_ph=DEFAULT_PH,
        #
        minimize=False,
        minimization_random_seed=DEFAULT_MINI_SEED,
        minimization_max_iterations=DEFAULT_MAX_IT,
        minimized_output=DEFAULT_MINI_OUT,
        #
        contact_cutoff=DEFAULT_CONTACT_CUTOFF,
        energy_log=DEFAULT_ENERGY_LOG,
        #
        silence_log=False,
        ):
    """
    Calculates interface atom-pairwise energies (LJ and Coulomb).

    Parameters
    ----------
    structure_fpath : str or Path
        The path to the input structure. The structure should be a
        complex between two or more proteins, that is, it should contain
        two or more chains.

    ff_name : str
        The name of the forcefield according to OpenMM.

    add_hydrogens : boolean
        Whether or not to add hydrogens according to the force field
        selected. Defaults to `False` (expects hydrogens to be already
        placed in the input structure).

    hydrogen_ph : float
        The pH value when adding hydrogens. Defaults to 7.0.

    hydrogen_variants : list
        A list of residue protonation states according to
        `openmm.app.modeller.Modeller.addHydrogens`. Defaults to
        `None`.

    minimize : boolean
        Whether to perform a short minimization in vaccum.

    minimization_random_seed : int
        Initial random seed number for the minimization. Defaults to
        `True`.

    minimization_max_iterations : int
        The number of iterations during the minimization according to
        `openmm.app.simulation.Simulation.minimizeEnergy`. Defaults to
        100.

    minimized_output : boolean or string or Path.
        `False` to not write the minimization output to a file.
        Otherwise, provide the name of the path. Writes the PDB file to
        the same folder as the `structure_fpath`.

    contact_cutoff : float
        The maximum distance to consider two atoms are in contact.
        Defaults to 5 Angstroms.

    energy_log : str, Path or `False`
        Name of the file to log the atom pairwise energies.
        To avoid writing the log file, provide `False`.

    silence_log : bool
        Whether to disable logging operations to the console. This may
        be useful if running this function in high throughput.

    Returns
    -------
    list of tuples
        Each item of the list is a 10 element tuple where the first 4
        indexes identify atom A the second for indexes identify atom B
        and the last two indexes are the LJ and Coulomb energies,
        respectively.

    Saves
    -----
        If `energy_log` is `True` saves to disk a text file with the
        atom pairs energies and the sum of all energies.
    """
    kcal_mole = units.kilocalorie_per_mole

    if silence_log:
        logging.disable(level=logging.INFO)
    logging.info('Started')

    logging.info('Reading structure: {}'.format(structure_fpath))

    # reads structure
    with open(structure_fpath, 'r') as fin:
        pdb = app.PDBFile(fin)

    # prepares forcefield
    forcefield = app.ForceField(ff_name)

    # makes a model of the structure
    model = app.Modeller(pdb.topology, pdb.positions)

    # Adds/replaces hydrogens according to user
    if add_hydrogens:
        logging.info('Adding hydrogens...')

        # first deletes all hydrogens
        to_delete = []
        for atom in pdb.topology.atoms():
            if atom.element == ElementHydrogen:
                to_delete.append(atom)
        model.delete(to_delete)
        del to_delete

        hydrogen_variants = _process_hydrogen_variants(hydrogen_variants)

        # add hydrogens
        model.addHydrogens(
            forcefield=forcefield,
            pH=hydrogen_ph,
            variants=hydrogen_variants,
            )

        # updates positions and topology (in case hydrogens were updated)
        pdb.topology = model.topology
        pdb.positions = model.positions

    logging.info('Preparing the system.')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)

    #
    # Setup simulation
    #
    integrator = mm.LangevinIntegrator(
        310 * units.kelvin,  # the actual value does not matter for the minimization
        1.0 / units.picosecond,
        2.0 * units.femtosecond,
        )
    integrator.setRandomNumberSeed(minimization_random_seed)
    integrator.setConstraintTolerance(0.00001)

    simulation = app.Simulation(pdb.topology, system, integrator)

    context = simulation.context
    context.setPositions(pdb.positions)

    #
    # Minimize
    #
    if minimize:
        logging.info('Minimizing potential energy.')

        state = context.getState(getEnergy=True)
        ini_ene = state.getPotentialEnergy().value_in_unit(kcal_mole)

        simulation.minimizeEnergy(maxIterations=minimization_max_iterations)
        pdb.positions = context.getState(getPositions=True).getPositions()

        state = context.getState(getEnergy=True)
        min_ene = state.getPotentialEnergy().value_in_unit(kcal_mole)
        logging.info('Pot. Energy: {:6.3f} (was {:6.3f}) kcal/mol'.format(min_ene, ini_ene))

        if minimized_output:
            _path = Path(structure_fpath)
            output_fname = Path(_path.parent, _path.stem + '_minimized.pdb')
            write_structure(output_fname, pdb)

    #
    # Store FF parameters for nonbonded force
    #
    # finds the nonbonded forces, those will be reset and used to calculate
    # the atom pairwise energies.
    logging.info('Resetting forces.')
    for f_idx in range(system.getNumForces()):
        f = system.getForce(f_idx)
        if isinstance(f, mm.NonbondedForce):
            nonb = f
            break

    atom_parameters = []
    for atom_idx in range(nonb.getNumParticles()):
        q, sigma, epsilon = nonb.getParticleParameters(atom_idx)
        atom_parameters.append((q, sigma, epsilon))

    del nonb

    #
    # Remove all forces acting on system
    #
    while system.getNumForces():
        system.removeForce(system.getNumForces() - 1)

    n_particles = system.getNumParticles()

    #
    # Configures LJ and Coulomb energies
    #
    # LJ 6-12
    logging.info('Setting up LJ(12-6) potential.')

    expr = "4*epsilon*((sigma/r)^12-(sigma/r)^6);"
    expr += "sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)"

    customLJForce = mm.CustomNonbondedForce(expr)
    customLJForce.addPerParticleParameter("sigma")
    customLJForce.addPerParticleParameter("epsilon")
    customLJForce.setNonbondedMethod(NONB_METHOD)
    customLJForce.setCutoffDistance(R_IJ)

    for particle in range(n_particles):
        _, sigma, epsilon = atom_parameters[particle]  # confirmation
        customLJForce.addParticle([_zero_sig, _zero_eps])

    # Coulomb
    logging.info('Setting up Coulomb potential.')

    ONE_4PI_EPS0 = 138.93563947857788
    expr = "({:g}*q)/r; q=q1*q2".format(ONE_4PI_EPS0)
    customCoulForce = mm.CustomNonbondedForce(expr)
    customCoulForce.addPerParticleParameter("q")
    customCoulForce.setNonbondedMethod(NONB_METHOD)
    customCoulForce.setCutoffDistance(R_IJ)

    for particle in range(n_particles):
        q, _, _ = atom_parameters[particle]
        customCoulForce.addParticle([_zero_q])

    #
    # Identifies atoms of each chain
    #
    atoms_in_chain = []
    for chain in pdb.topology.chains():
        atoms_in_chain.append([])
        for atom in chain.atoms():
            atoms_in_chain[-1].append(atom.index)
    _chains_lens = [len(x) for x in atoms_in_chain]
    logging.info(f'Lengths of chains: {_chains_lens}.')

    system.addForce(customLJForce)
    system.addForce(customCoulForce)

    lj = system.getForce(0)
    coul = system.getForce(1)

    lj.setForceGroup(1)
    coul.setForceGroup(2)

    context.reinitialize(preserveState=True)

    # Atom-Atom
    all_atoms = list(pdb.topology.atoms())
    n_atoms = len(all_atoms)
    logging.info(f'System has a total of {n_atoms} atoms.')

    # OpenMM numbers residues consecutively regardless of the chain ID.
    # Here we calculate the number of residues per chain such that on
    # the output the residues are numbered from 1 starting on each
    # chain. In the resulting dictionary keys are chain IDs and values
    # the number to subtract to the OpenMM residue numeration.
    residues_per_chain = {}
    for chain in pdb.topology.chains():
        residues_per_chain[chain] = len(list(chain.residues()))
    residues_offset = list(
        itertools.accumulate(
            list(residues_per_chain.values())[:-1],
            operator.add)
            )
    residues_offset = [0] + residues_offset
    for key, value in zip(residues_per_chain.keys(), residues_offset):
        residues_per_chain[key] = value - 1

    #
    # Calculate atom pairwise energy for the complex interface
    #
    # Usually the input PDB will have two chains. But the algorithm is
    # ready to calculate pairwise combination for multimeric complexes.
    logging.info('Calculating energies...')

    lj_tot, coul_tot, ener_log = 0.0, 0.0, []
    sq_min_cutoff = contact_cutoff * contact_cutoff

    for chain_A, chain_B in itertools.combinations(list(pdb.topology.chains()), r=2):
        _atoms_of_both_chains = list(chain_A.atoms()) + list(chain_B.atoms())

        for atom_i, atom_j in itertools.combinations(_atoms_of_both_chains, r=2):

            # skips hydrogens and atoms of the same chain
            if atom_i.residue.chain == atom_j.residue.chain:
                continue

            atom_i_idx = atom_i.index
            atom_j_idx = atom_j.index

            xyz_i = pdb.positions[atom_i_idx]
            xyz_j = pdb.positions[atom_j_idx]
            distance_square = sq_eucl_dij(xyz_i, xyz_j)
            if distance_square > sq_min_cutoff:  # contact
                continue

            q_i, sig_i, eps_i = atom_parameters[atom_i_idx]

            # set parameters for atom i
            lj.setParticleParameters(atom_i_idx, [sig_i, eps_i])
            lj.updateParametersInContext(context)

            coul.setParticleParameters(atom_i_idx, [q_i])
            coul.updateParametersInContext(context)

            # set parameters for atom j
            q_j, sig_j, eps_j = atom_parameters[atom_j_idx]
            lj.setParticleParameters(atom_j_idx, [sig_j, eps_j])
            lj.updateParametersInContext(context)

            coul.setParticleParameters(atom_j_idx, [q_j])
            coul.updateParametersInContext(context)

            # Calculate energies
            state = context.getState(getEnergy=True, groups={1})
            lj_ene = state.getPotentialEnergy().value_in_unit(kcal_mole)

            state = context.getState(getEnergy=True, groups={2})
            coul_ene = state.getPotentialEnergy().value_in_unit(kcal_mole)

            # 'Turn off' first particle
            lj.setParticleParameters(atom_i_idx, [_zero_sig, _zero_eps])
            lj.updateParametersInContext(context)

            coul.setParticleParameters(atom_i_idx, [_zero_q])
            coul.updateParametersInContext(context)

            # 'Turn off' second particle
            lj.setParticleParameters(atom_j_idx, [_zero_sig, _zero_eps])
            lj.updateParametersInContext(context)

            coul.setParticleParameters(atom_j_idx, [_zero_q])
            coul.updateParametersInContext(context)

            if lj_ene != 0.0 or coul_ene != 0.0:
                lj_tot += lj_ene
                coul_tot += coul_ene

                pair_report = (
                        atom_i.residue.chain.id,
                        atom_i.residue.name,
                        int(atom_i.residue.index) - residues_per_chain[atom_i.residue.chain],
                        atom_i.name,
                        atom_j.residue.chain.id,
                        atom_j.residue.name,
                        int(atom_j.residue.index) - residues_per_chain[atom_j.residue.chain],
                        atom_j.name,
                        lj_ene,
                        coul_ene,
                        )

                ener_log.append(pair_report)

    if energy_log:
        logging.info(f'Saving energy log to {energy_log!r}.')
        # chain, resname, res number, atom name
        log_fmt = '{:<2}{:<4}{:<5}{:<4}'
        log_fmt = log_fmt + ' - ' + log_fmt + ' {:+11.5f} {:+11.5f}'
        _ = ['chainA resnameA resiA atomA - chainB resnameB resiB atomA LJ Coulomb (kcal/mole)']
        _.extend((log_fmt.format(*pair) for pair in ener_log))
        _.append(f'Total LJ: {lj_tot:.5f} (kcal/mole)')
        _.append(f'Total Coulomb: {coul_tot:.5f} (kcal/mole)')

        with open(energy_log, 'w') as fout:
            fout.writelines(os.linesep.join(_))

    return ener_log


def load_args(ap):
    cmd = ap.parse_args()
    return cmd


def main_cli():
    cmd = load_args(ap)
    main(**vars(cmd))


def main(*args, **kwargs):
    calculate_interface_atom_pairwise_energies(*args, **kwargs)


if __name__ == '__main__':
    main_cli()
