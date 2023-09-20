#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014-2023 pdb-tools project, https://github.com/haddocking/pdb-tools
# Copyright 2023 João M.C. Teixeira (@joaomcteixeira).
#
# https://www.bonvinlab.org/pdb-tools/
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# For the `pdbprep` project, pdb-tools scripts were modified to reduce their
# versatility in exchange for additonal speed. The core algorithm of the script
# was maintained, but boilerplate to create a versatile command-line was
# removed, and some arguments were blocked to specific options.  Modifications
# to the original `pdb-tools` scripts were done by João M.C. Teixeira
# (@joaomcteixeira) who is co-author in the `pdb-tools` project. Any credit to
# the value of `pdb-tools` scripts present in the `pdbprep` project should be
# given to the original `pdb-tools` project.
import os
import sys


def run(fhandle, sorting_keys):
    # Sort keys
    chain_key = lambda x: x[21]  # chain id
    resid_key = lambda x: (int(x[22:26]), x[26])  # resid, icode
    atoms_key = lambda x: int(x[6:11])  # atom serial
    altloc_key = lambda x: x[16]
    icode_key = lambda x: x[26]

    # ignored fields
    ignored = (('END', 'MASTER', 'TER'))

    # First, separate records
    header_data = []
    atomic_data = []
    hetatm_data = []
    anisou_data = {}  # Matches a unique atom uid
    conect_data = []
    chains_list = []  # list of chain identifiers in order they were read
    for line in fhandle:
        if line.startswith('ATOM'):
            atomic_data.append(line)
            chains_list.append(line[21])
        elif line.startswith('HETATM'):
            hetatm_data.append(line)
            chains_list.append(line[21])
        elif line.startswith(ignored):
            continue
        elif line.strip():  # remove empty lines
            header_data.append(line)

    # Map original chain orders to integers
    # To circumvent mixed chains when sorting by residue number
    chain_order = {}
    ordinal = 0
    for chain in chains_list:
        if chain not in chain_order:
            chain_order[chain] = ordinal
            ordinal += 1
    del chains_list

    # Always sort by altloc
    atomic_data.sort(key=atoms_key)
    atomic_data.sort(key=altloc_key)
    hetatm_data.sort(key=atoms_key)
    hetatm_data.sort(key=altloc_key)

    if 'R' in sorting_keys:
        atomic_data.sort(key=icode_key)
        atomic_data.sort(key=resid_key)
        hetatm_data.sort(key=icode_key)
        hetatm_data.sort(key=resid_key)

    if 'C' in sorting_keys:
        atomic_data.sort(key=chain_key)
        hetatm_data.sort(key=chain_key)
    else:  # restore original order
        chain_key = lambda x: chain_order.get(x[21])
        atomic_data.sort(key=chain_key)
        hetatm_data.sort(key=chain_key)

    # Sort conect statements by the central atom
    # Share the same format at ATOM serial number
    conect_data.sort(key=atoms_key)

    # Now return everything in order:
    #  - ATOMs intercalated with ANISOU
    #  - HETATM
    #  - CONECT
    sorted_data = header_data + atomic_data + hetatm_data + conect_data
    for line in sorted_data:

        yield line

        atom_uid = line[12:27]
        anisou_record = anisou_data.get(atom_uid)
        if anisou_record:
            yield anisou_record



def main():
    # Check Input
    pdbfh = sys.stdin

    # Do the job
    new_pdb = run(pdbfh, 'CR')

    try:
        _buffer = []
        _buffer_size = 5000  # write N lines at a time
        for lineno, line in enumerate(new_pdb):
            if not (lineno % _buffer_size):
                sys.stdout.write(''.join(_buffer))
                _buffer = []
            _buffer.append(line)

        sys.stdout.write(''.join(_buffer))
        sys.stdout.flush()
    except IOError:
        # This is here to catch Broken Pipes
        # for example to use 'head' or 'tail' without
        # the error message showing up
        pass

    # last line of the script
    # We can close it even if it is sys.stdin
    pdbfh.close()
    sys.exit(0)


if __name__ == '__main__':
    main()
