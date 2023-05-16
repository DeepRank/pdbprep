#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys


def run(fhandle, resname_set):
    """
    Remove specific residue that match a given name.

    Non-coords lines are maintained.

    This function is a generator.

    Parameters
    ----------
    fhandle : a line-by-line iterator of the original PDB file.

    Yields
    ------
    str (line-by-line)
        The PDB lines not matching the residues selected.
        Non-coord lines are yielded as well.
    """
    records = ('ATOM', 'HETATM', 'ANISOU', 'TER')
    for line in fhandle:
        if line.startswith(records):
            if line[17:20].strip() in resname_set:
                continue
        yield line


def main():
    # Check Input
    pdbfh, resname_set = sys.stdin, ('HOH',)

    # Do the job
    new_pdb = run(pdbfh, resname_set)

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
