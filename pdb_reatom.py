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


def run(fhandle, starting_value):
    # CONECT 1179  746 1184 1195 1203
    fmt_CONECT = "CONECT{:>5s}{:>5s}{:>5s}{:>5s}{:>5s}" + " " * 49 + os.linesep
    char_ranges = (slice(6, 11), slice(11, 16),
                   slice(16, 21), slice(21, 26), slice(26, 31))

    serial_equiv = {'': ''}  # store for conect statements

    serial = starting_value
    records = ('ATOM', 'HETATM')
    for line in fhandle:
        if line.startswith(records):
            serial_equiv[line[6:11].strip()] = serial
            yield line[:6] + str(serial).rjust(5) + line[11:]
            serial += 1
            if serial > 99999:
                emsg = 'Cannot set atom serial number above 99999.\n'
                sys.stderr.write(emsg)
                sys.exit(1)



def main():

    # Do the job
    new_pdb = run(sys.stdin, 1)

    # Output results
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
    # Close file handle even if it is sys.stdin, no problem here.
    #pdbfh.close()
    sys.exit(0)


if __name__ == '__main__':
    main()
