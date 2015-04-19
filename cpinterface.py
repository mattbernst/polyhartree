# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import sys
from sharedutilities import Utility
import geoprep

try:
    from ebsel import EMSL_api
except ImportError:
    msg = "Warning! Unable to import ebsel. Attempts to prepare calculations that need basis set data will fail. Make sure ebsel is on your PYTHONPATH to run calculations using basis sets. https://github.com/mattbernst/ebsel\n"
    sys.stderr.write(msg)

class Messages(object):
    def log(self, msg):
        self.messages.append(msg)

    def log_once(self, msg):
        if not self.messages or self.messages[-1] != msg:
            self.messages.append(msg)

class Job(Utility, Messages):
    def __init__(self, deck="", system=None, runstate="begin", extras={}):
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
        self.extras = extras

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

    def make_energy_job(self, system, method, options={}):
        raise NotImplementedError

    def fragment_to_system(self, sf):
        """It's convenient to be able to pass in a single-fragment system as
        a fragment without manually converting it to a system.

        @param sf: a system or a single fragment
        @type sf : geoprep.Fragment | geoprep.System
        @return: a system
        @rtype : geoprep.System
        """

        if type(sf) != geoprep.System:
            s = geoprep.System(sf)
            return s
        else:
            return sf

    def make_opt_job(self, system, method, options={}):
        raise NotImplementedError

    def get_basis_data(self, system, options={}):
        """Get basis set data for all atoms, and raise errors for missing
        or incorrectly specified basis sets.

        options:
         basis_format: "gamess-us", "nwchem", or "gaussian94"

        @param system: molecular system
        @type system : geoprep.System
        @param options: ignored
        @type options : dict
        @return: basis data block
        @rtype : str
        """

        basis_format = options["basis_format"]
        basis_names = system.atom_properties("basis_name")
        symbols = system.atom_properties("symbols")
        groups = {}
        basis_groups = {}
        function_types = {}

        if None in basis_names:
            missing = []
            for j, name in enumerate(basis_names):
                if name is None:
                    missing.append(j)
                    
            raise ValueError("Missing basis set names on atoms {0}".format(missing))

        for j, name in enumerate(basis_names):
            symbol = symbols[j]
            try:
                groups[name].add(symbol)
            except KeyError:
                groups[name] = set([symbol])

        el = EMSL_api.get_EMSL_local(fmt=basis_format)
        for group in groups:
            provided_elements = set(el.get_list_element_available(group))

            #no provided elements at all means unknown basis set name
            if not provided_elements:
                raise ValueError("Unknown basis set {0}".format(group))
                
            elements = groups[group]
            missing_elements = elements.difference(provided_elements)

            if missing_elements:
                raise ValueError("Basis set {0} missing parameters for elements {1}".format(repr(group), list(missing_elements)))

            basis_data = {}
            for element in elements:
                bd = el.get_basis(group, [element])[0]
                basis_data[element] = bd
                
            basis_groups[group] = basis_data

            soc = el.spherical_or_cartesian(group)
            try:
                function_types[soc].append(group)
            except KeyError:
                function_types[soc] = [group]

        #can't mix spherical and cartesian basis sets in a single system
        if len(function_types) > 1:
            raise ValueError("Attempted mixing spherical and cartesian basis sets: {0}".format(function_types))

        #note the basis function type as an annotation
        d = {"spherical_or_cartesian" : function_types.keys()[0],
             "data" : basis_groups}
        
        return d

    def check_element_support(self, system, options={}):
        raise NotImplementedError

    def check_electronic_reference(self, reference):
        """Check that electronic reference is supported for calculation to
        be attempted. Raise an exception on unsupported reference.

        Values that might be supported by individual packages/methods include
        RHF, UHF, ROHF for wavefunction methods and RDFT, UDFT, RODFT, UKS,
        ROKS, RKS for density functional methods.

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
        """Check that calculation method is supported by the computational
        chemistry back-end in use.

        Some examples: semiempirical:pm3, semiempirical:am1, hf:uhf,
        dft:m05:roks, correlated:mp2:rhf, mm:amber-99...

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
