#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys


def run(fhandle, strict=False):
    not_strict = not strict
    fhandle = iter(fhandle)

    def make_TER(prev_line):
        """Creates a TER statement based on the last ATOM/HETATM line.
        """

        # Add last TER statement
        serial = int(prev_line[6:11]) + 1
        rname = prev_line[17:20]
        chain = prev_line[21]
        resid = prev_line[22:26]
        icode = prev_line[26]

        return fmt_TER.format(serial, rname, chain, resid, icode)

    # TER     606      LEU A  75
    fmt_TER = "TER   {:>5d}      {:3s} {:1s}{:>4s}{:1s}" + " " * 53 + "\n"

    records = ('ATOM', 'HETATM')
    ignored = ('TER', 'END', 'CONECT', 'MASTER', 'ENDMDL')
    # Iterate up to the first ATOM/HETATM line
    prev_line = None
    num_models = 1
    in_model = False
    for line in fhandle:

        line = line.strip()  # We will pad/add \n later to make uniform

        if line.startswith('MODEL'):
            line = "MODEL " + "    " + str(num_models).rjust(4)
            num_models += 1
            in_model = True

        if line.startswith(ignored):  # to avoid matching END _and_ ENDMDL
            continue

        # Check line length
        line = "{:<80}\n".format(line)

        yield line

        if line.startswith(records):
            prev_line = line
            break

    # Now go through all the remaining lines
    atom_section = False
    serial_offset = 0  # To offset after adding TER records
    for line in fhandle:

        line = line.strip()

        if line.startswith(ignored):
            continue

        # Treat ATOM/HETATM differently
        #   - no TER in HETATM
        if line.startswith('ATOM'):

            is_gap = (int(line[22:26]) - int(prev_line[22:26])) > 1
            if atom_section and (line[21] != prev_line[21] or (not_strict and is_gap)):
                serial_offset += 1  # account for TER statement
                yield make_TER(prev_line)

            serial = int(line[6:11]) + serial_offset
            line = line[:6] + str(serial).rjust(5) + line[11:]
            prev_line = line
            atom_section = True

        elif line.startswith('HETATM'):
            if atom_section:
                atom_section = False
                serial_offset += 1  # account for TER statement
                yield make_TER(prev_line)

            serial = int(line[6:11]) + serial_offset
            line = line[:6] + str(serial).rjust(5) + line[11:]
            prev_line = line

        elif line.startswith('ANISOU'):
            # Fix serial based on previous atom
            # Avoids doing the offset again
            serial = int(prev_line[6:11])
            line = line[:6] + str(serial).rjust(5) + line[11:]

        else:
            if atom_section:
                atom_section = False
                yield make_TER(prev_line)
                if in_model:
                    yield "{:<80}\n".format("ENDMDL")
                    in_model = False

            if line.startswith('MODEL'):
                line = "MODEL " + "    " + str(num_models).rjust(4)
                num_models += 1
                in_model = True
                serial_offset = 0

        if serial > 99999:
            emsg = 'ERROR!! Structure contains more than 99.999 atoms.\n'
            sys.stderr.write(emsg)
            sys.stderr.write(__doc__)
            sys.exit(1)

        # Check line length
        line = "{:<80}\n".format(line)

        yield line

    else:
        if atom_section:
            # Add last TER statement
            atom_section = False
            yield make_TER(prev_line)
            if in_model:
                yield "{:<80}\n".format("ENDMDL")
                in_model = False

    # Add END statement
    yield "{:<80}\n".format("END")



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
