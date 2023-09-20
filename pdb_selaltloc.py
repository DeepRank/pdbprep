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


def run(fhandle):
    records = ('ATOM', 'HETATM')
    #terminators = ('TER', 'END', 'CONECT', 'END', 'ENDMDL', 'MODEL')
    #meaningful = records + terminators

    # register atom information
    register = dict()

    # register comment lines
    others = []

    # register current chain
    chain = None
    prev_chain = None

    # keep record of the line number. This will be used to sort lines
    # after selecting the desired alternative location
    nline = 0

    # the loop will collect information on the different atoms
    # throughout the PDB file until a new chain or any terminal line is
    # found. At that point, the collected information is flushed because
    # all altlocs for that block have been defined.
    for line in fhandle:
        nline += 1

        if line.startswith(records):

            # here resnum + insertion code are taken to identify
            # different residues
            resnum = line[22:27]
            atomname = line[12:16]
            altloc = line[16]
            chain = line[21:22]

            # flush lines because we enter a new chain
            if chain != prev_chain:
                # the "yield from" statement is avoided to keep
                # compatibility with Python 2.7
                for _line in _flush(register):
                    yield _line

                # Python 2.7 compatibility. Do not use .clear() method
                # restart help variables
                del register, others
                register, others = dict(), []

            # organizes information hierarchically
            resnum_d = register.setdefault(resnum, {})
            atomname_d = resnum_d.setdefault(atomname, {})
            altloc_d = atomname_d.setdefault(altloc, [])

            # adds info to dictionary
            altloc_d.append((nline, line))

        # flush information because we reached the end of a block
        #elif line.startswith(terminators):
        #    for _line in _flush(register):
        #        yield _line

        #    del register, others
        #    register, others = dict(), []

        #    yield line  # yield the current line after flush

        prev_chain = chain

    # at the end of the PDB, flush the remaining lines
    for _line in _flush(register):
        yield _line


def _flush(register):
    """
    Processes the collected atoms according to the selaltloc option.
    """
    lines_to_yield = []

    atom_lines = ('ATOM', 'HETATM')

    # anisou lines are treated specially
    anisou_lines = ('ANISOU',)

    for resnum, atomnames in register.items():

        for atomname, altlocs in atomnames.items():

            # gathers all alternative locations for the atom
            all_lines = []
            for altloc, lines in altlocs.items():
                all_lines.extend(lines)

            # identifies the highest occupancy combining dictionary
            # and sorting
            new = {}
            for line_number, line in all_lines:
                if line.startswith(atom_lines):
                    occupancy_number = line[54:60]
                    list_ = new.setdefault(occupancy_number, [])
                    list_.append((line_number, line))

                # assumes ANISOU succeed the respective ATOM line
                elif line.startswith(anisou_lines):
                    list_.append((line_number, line))

            # sort keys by occupancy
            keys_ = sorted(new.keys(), key=lambda x: float(x.strip()), reverse=True)

            these_atom_lines = new[keys_[0]]
            if len(keys_) == 1 and len(these_atom_lines) > 1:
                # address "take first if occ is the same"
                # see: https://github.com/haddocking/pdb-tools/issues/153#issuecomment-1488627668
                lines_to_yield.extend(_remove_altloc(these_atom_lines[0:1]))

                # if there's ANISOU, add it
                if these_atom_lines[1][1].startswith(anisou_lines):
                    lines_to_yield.extend(_remove_altloc(these_atom_lines[1:2]))

            # this should run when there are more than one key or
            # the key has only one atom line. Keys are the occ
            # value.
            else:
                # when occs are different, select the highest one
                lines_to_yield.extend(_remove_altloc(these_atom_lines))

            del all_lines, new


    # lines are sorted to the line number so that the output is sorted
    # the same way as in the input PDB
    lines_to_yield.sort(key=lambda x: x[0])

    # the line number is ignored, only the line is yield
    for line_number, line in lines_to_yield:
        yield line


def _remove_altloc(lines):
    # the altloc ID is removed in processed altloc lines
    for line_num, line in lines:
        yield (line_num, line[:16] + ' ' + line[17:])


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
