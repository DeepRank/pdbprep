#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys


def pad_line(line):
    """Helper function to pad line to 80 characters in case it is shorter"""
    size_of_line = len(line)
    if size_of_line < 80:
        padding = 80 - size_of_line + 1
        line = line.strip('\n') + ' ' * padding + '\n'
    return line[:81]  # 80 + newline character


def run(fhandle, starting_resid):
    """
    Reset the residue number column to start from a specific number.

    This function is a generator.

    Parameters
    ----------
    fhandle : a line-by-line iterator of the original PDB file.

    starting_resid : int
        The starting residue number.

    Yields
    ------
    str (line-by-line)
        The modified (or not) PDB line.
    """
    _pad_line = pad_line
    prev_resid = None  # tracks chain and resid
    resid = starting_resid - 1  # account for first residue
    records = ('ATOM', 'HETATM')
    for line in fhandle:
        line = _pad_line(line)
        if line.startswith(records):
            line_resuid = line[17:27]
            if line_resuid != prev_resid:
                prev_resid = line_resuid
                resid += 1
                if resid > 9999:
                    emsg = 'Cannot set residue number above 9999.\n'
                    sys.stderr.write(emsg)
                    sys.exit(1)

            yield line[:22] + str(resid).rjust(4) + line[26:]



def main():
    # Check Input

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
