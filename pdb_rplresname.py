#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys


def run(fhandle):
    replacements = {
        "MSE": "MET",
        "HIP": "HIS",
        "HIE": "HIS",
        "HID": "HIS",
        "HSE": "HIS",
        "HSD": "HIS",
        }
    for line in fhandle:
        line_resname = line[17:20].strip()
        name_to = replacements.get(line_resname, line_resname)
        yield line[:17] + name_to.rjust(3) + line[20:]



def main():
    # Check Input
    pdbfh = sys.stdin

    # Do the job
    new_pdb = run(pdbfh)

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
    pdbfh.close()
    sys.exit(0)


if __name__ == '__main__':
    main()
