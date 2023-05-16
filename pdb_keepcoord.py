#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
