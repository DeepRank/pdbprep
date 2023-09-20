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
    records = ('ATOM  ', 'HETATM')
    for line in fhandle:
        if line.startswith(records):
            yield line


def main():
    # Check Input
    pdbfh = open(sys.argv[1], 'r')

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
