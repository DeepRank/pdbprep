#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014-2023 pdb-tools project.
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


def run(fhandle):
    option_set = set()  # empty if option_list is empty

    # Keep track of residue numbering
    # Keep track of residues read (chain, resname, resid)
    offset = 0
    prev_resi = None
    seen_ids = set()
    clean_icode = False
    records = ('ATOM', 'HETATM', 'ANISOU', 'TER')
    for line in fhandle:

        if line.startswith(records):
            res_uid = line[17:27]  # resname, chain, resid, icode
            id_res = line[21] + line[22:26].strip()  # A99, B12
            has_icode = line[26].strip()  # ignore ' ' here

            # unfortunately, this is messy but not all PDB files follow a nice
            # order of ' ', 'A', 'B', ... when it comes to insertion codes..
            if prev_resi != res_uid:  # new residue
                # Does it have an insertion code
                # OR have we seen this chain + resid combination before?
                # #2 to catch insertions WITHOUT icode ('A' ... ' ' ... 'B')
                if (has_icode or id_res in seen_ids):
                    # Do something about it
                    # if the user provided options and this residue is in them
                    # OR if the user did not provide options
                    if (not option_set) or (id_res in option_set):
                        clean_icode = True
                    else:
                        clean_icode = False
                else:
                    clean_icode = False

                prev_resi = res_uid

                if id_res in seen_ids:  # offset only if we have seen this res.
                    offset += 1

            if clean_icode:  # remove icode
                line = line[:26] + ' ' + line[27:]

            # Modify resid if necessary
            resid = int(line[22:26]) + offset
            line = line[:22] + str(resid).rjust(4) + line[26:]
            seen_ids.add(id_res)

            # Reset offset on TER
            if line.startswith('TER'):
                offset = 0

        yield line



def main():
    # Check Input
    pdbfh = sys.stdin

    # Do the job
    new_pdb = run(pdbfh)

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
