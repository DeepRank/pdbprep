#!/usr/bin/env python
import sys
import json


def detect_protonation_state(resname, residues, atoms_in_residue):
    if resname == 'HIS':
        if 'HD1' in atoms_in_residue and 'HE2' in atoms_in_residue:
            resname = 'HIP'
        elif 'HD1' in atoms_in_residue:
            resname = 'HID'
        elif 'HE2' in atoms_in_residue:
            resname = 'HIE'
        else:
            resname = 'HIN'
    elif resname == 'ASP':
        if 'HD2' in atoms_in_residue or 'HD1' in atoms_in_residue:
            resname = 'ASN'
    elif resname == 'GLU':
        if 'HE2' in atoms_in_residue or 'HE1' in atoms_in_residue:
            resname = 'GLH'
    elif resname == 'LYS':
        if not all(_a in atoms_in_residue for _a in ('HZ1', 'HZ2', 'HZ3')):
            resname = 'LYN'
    elif resname == 'CYS':
        if 'HG' in atoms_in_residue:
            resname = 'CYS'
        else:
            resname = 'CYX'
    else:
        raise AssertionError(
            f'Found unexpected residue: {resname}. '
            'Code shouldn\'t be here'
            )

    residues.append(resname)
    return


fin_name = sys.argv[1]
fout_name = sys.argv[2]
fin = open(fin_name, 'r')

records = ('ATOM',)
protonable_residues = ('HIS', 'ASP', 'GLU', 'CYS', 'LYS')
prev_resid = None
atoms_in_residue = set()
residues = []
resname = None
for line in fin:

    if not line.startswith(records):
        continue

    resid = line[21:26]  # chain ID + res number
    resname = line[17:20]

    if resid != prev_resid and len(atoms_in_residue) > 4:

        if prev_resname in protonable_residues:
            detect_protonation_state(prev_resname, residues, atoms_in_residue)

        else:
            residues.append(None)

        atoms_in_residue.clear()

    atom_name = line[12:16].strip()
    atoms_in_residue.add(atom_name)

    prev_resid = resid
    prev_resname = resname

else:
    if resname in protonable_residues:
        detect_protonation_state(prev_resname, residues, atoms_in_residue)
    else:
        residues.append(None)

fin.close()

with open(fout_name, 'w') as fout:
    json.dump({fin_name: residues}, fout, indent=4)
