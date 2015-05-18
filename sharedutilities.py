# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import subprocess
import shlex

ELEMENTS = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
            "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
            "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
            "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
            "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
            "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
            "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
            "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
            "Md", "No", "Lr"]

class Utility(object):
    def execute_local(self, cmd, stdin_data="", bash_shell=False, cwd=None):
        """Execute a command with subprocess.Popen, optionally supplying
        data to the command through stdin, and return the results.

        :param cmd: a command line for an external program
        @type cmd : str
        :param stdin_data: optional data to supply on stdin to external program
        @type stdin_data : str
        :param bash_shell: execute command through Bash with rcfile, if True
        @type bash_shell : bool
        @return: (data from stdout, return code)
        @rtype : tuple
        """

        command = shlex.split(cmd)
        with(open("/dev/null", "w")) as devnull:
            if bash_shell:
                command = ["/bin/bash", "-i", "-c"] + [" ".join(command)]
                
            p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                 stdin=subprocess.PIPE, stderr=devnull,
                                 cwd=cwd)
            output = p.communicate(input=stdin_data)[0]

        return (output, p.returncode)

    def numericize(self, line, numeric_only=False):
        """Split a line of text by whitespace and try to convert all
        numbers to floating point numeric values.

        e.g.

        "Root 1 singlet 19.48009 a.u. 530.08030 eV"

        becomes

        ["Root", 1.0, "singlet", 19.48009, "a.u.", 530.0803, "eV"]

        :param line: input line of text from a nwparse QA file
        @type line : str
        :param numeric_only: exclude unconverted fragments if True
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

    def n_number_from_line(self, line, n, m):
        """Get the Nth numeric value from a line where M values are present.

        :param line: a line from a data file
        @type line : str
        @return: found numeric value
        @rtype : float
        """

        numbers = self.numericize(line, numeric_only=True)
        if len(numbers) != m:
            raise ValueError("Expected {0} numeric values in line: {1}".format(m, line))
        return numbers[n]

    def ev_to_au(self, v):
        """Convert an energy expressed in electron volts to a Hartree value.

        :param v: energy in electron volts
        @type v : float
        @return: energy in Hartrees
        @rtype : float
        """

        return v / 27.211385

    def kcalm_to_au(self, v):
        """Convert an energy expressed in kcal/mol to a Hartree value.

        :param v: energy in kcal/mol
        @type v : float
        @return: energy in Hartrees
        @rtype : float
        """

        return v / 627.509469

    def bohr_to_angstrom(self, v):
        """Convert a distance expressed in bohr to angstrom.

        :param v: distance in bohr
        @type v : float
        @return: distance in angstrom
        @rtype : float
        """

        return v * 0.52917721092
