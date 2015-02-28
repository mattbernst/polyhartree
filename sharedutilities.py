# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import subprocess
import shlex

class Utility(object):
    def execute(self, cmd, stdin_data='', shell=False):
        """Execute a command with subprocess.Popen, optionally supplying
        data to the command through stdin, and return the results.

        @param cmd: a command line for an external program
        @type cmd : str
        @param stdin_data: optional data to supply on stdin to external program
        @type stdin_data : str
        @param shell: execute command through system shell, if True
        @type shell : bool
        @return: (data from stdout, return code)
        @rtype : tuple
        """

        command = shlex.split(cmd)
        with(open('/dev/null', 'w')) as devnull:
            p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                 stdin=subprocess.PIPE, stderr=devnull,
                                 shell=shell)
            output = p.communicate(input=stdin_data)[0]

        return (output, p.returncode)

    def numericize(self, line, numeric_only=False):
        """Split a line of text by whitespace and try to convert all
        numbers to floating point numeric values.

        e.g.

        'Root 1 singlet 19.48009 a.u. 530.08030 eV'

        becomes

        ['Root', 1.0, 'singlet', 19.48009, 'a.u.', 530.0803, 'eV']

        @param line: input line of text from a nwparse QA file
        @type line : str
        @param numeric_only: exclude unconverted fragments if True
        @type numeric_only : bool
        @return: mixed list of strings and floats
        @rtype : list
        """

        converted = []
        for piece in line.split():
            try:
                v = float(piece)
            except ValueError:
                v = piece

            if type(v) == float or numeric_only == False:
                converted.append(v)

        return converted

    def one_number_from_line(self, line):
        """Get a single numeric value from a line where only one is present.

        @param line: a line from a data file
        @type line : str
        @return: found numeric value
        @rtype : float
        """

        numbers = self.numericize(line, numeric_only=True)
        if len(numbers) != 1:
            raise ValueError("Expected exactly one value in line: {0}".format(line))
        return numbers[0]

    def ev_to_au(self, v):
        """Convert an energy expressed in electron volts to a Hartree value.

        @param v: energy in electron volts
        @type v : float
        @return: energy in Hartrees
        @rtype : float
        """

        return v / 27.211385

    def kcalm_to_au(self, v):
        """Convert an energy expressed in kcal/mol to a Hartree value.

        @param v: energy in kcal/mol
        @type v : float
        @return: energy in Hartrees
        @rtype : float
        """

        return v / 627.509469 
