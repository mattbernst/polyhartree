# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-

import string
from sharedutilities import Utility

class Messages(object):
    def log(self, msg):
        self.messages.append(msg)

    def log_once(self, msg):
        if not self.messages or self.messages[-1] != msg:
            self.messages.append(msg)

class Job(Utility, Messages):
    def __init__(self, deck="", system=None, runstate="begin"):
        #states: begin, running, complete, error
        self.runstate = runstate
        self.system = system
        self.deck = deck
        self.stdout = ""
        self.logdata = ""
        self.energy = None
        self.heat_of_formation = None
        self.geometry = None
        self.messages = []

    def run_local(self, options={}):
        raise NotImplementedError

    def extract_last_energy(self, data, options={}):
        raise NotImplementedError

    def extract_heat_of_formation(self, data, options={}):
        raise NotImplementedError

class MolecularCalculator(Messages):
    def __init__(self, *args, **kw):
        self.messages = []
        
    def create_geometry(self, molecule, options={}):
        raise NotImplementedError

    def extract_geometry(self, data, options={}):
        raise NotImplementedError

    def extract_energy(self, data, options={}):
        raise NotImplementedError

    def make_energy_job(self, molecule, method, options={}):
        raise NotImplementedError

    def make_opt_job(self, molecule, method, options={}):
        raise NotImplementedError

    def check_element_support(self, system, options={}):
        raise NotImplementedError

    def get_elements(self, system):
        """Get all elements incorporated in system.

        @param system: molecular system
        @type system : cinfony molecule
        @return: unique list of system elements, in order of atomic number
        @rtype : list
        """
        
        elements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
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
        
        formula = system.formula
        for char in formula:
            if char not in string.ascii_letters:
                formula = formula.replace(char, " ")

        formula_elements = set(formula.split())
        included_elements = [e for e in elements if e in formula_elements]
        
        return included_elements

    def check_electronic_reference(self, reference):
        """Check that electronic reference is supported for calculation to
        be attempted. Raise an exception on unsupported reference.

        Values that might be supported by individual packages/methods include
        RHF, UHF, ROHF for wavefunction methods RDFT, UDFT, RODFT, UKS, ROKS,
        RKS for density functional methods.

        May need to do special things here for multi-reference methods.

        @param reference: reference name to check
        @type reference : str
        @return: original reference name
        @rtype : str
        """
        
        if reference in self.references:
            return reference

        else:
            raise ValueError("Unrecognized reference scheme {0}".format(repr(reference)))

    def check_method(self, method):
        """Check that calculation method is supported by the quantum chemistry
        back-end in use.

        Some examples: semiempirical:pm3, semiempirical:am1, hf:uhf, mp2:rhf...

        @param method: method name to check
        @type reference : str
        @return: original method name
        @rtype : str
        """
        
        if method in self.methods:
            return method
            
        else:
            raise ValueError("Unrecognized method {0}".format(repr(method)))

    def check_coordinates(self, coordinate_choice):
        """Check that coordinate system is supported by the quantum chemistry
        back-end in use.

        Some examples: cartesian, zmatrix

        @param coordinate_choice: coordinate system name to check
        @type reference : str
        @return: original coordinate system name
        @rtype : str
        """
        
        if coordinate_choice in self.coordinate_choices:
            return coordinate_choice

        else:
            raise ValueError("Unrecognized coordinates option {0}".format(repr(coordinate_choice)))
