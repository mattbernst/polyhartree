# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import hashlib
import uuid

import cpinterface
import geoprep
import sharedutilities

class Mopac7Job(cpinterface.Job):
    def __init__(self, *args, **kw):
        super(Mopac7Job, self).__init__(*args, **kw)
        self.backend = "mopac7"
        #intended to extract XYZ geometric coordinates from a line like
        #2         C                  1.3220     .0000     .0000
        self.geometry_matcher = ([int, str, float, float, float],
                                 [2, 3, 4])

    def extract_last_energy(self, data, options={}):
        """Get last energy message from log file and store it as self.energy.

        :param data: log file contents
        :type data : str
        :param options: ignored
        :type options : dict
        """

        electronic_energy = None
        core_repulsion = None
        
        for line in data.split("\n"):
            if "ELECTRONIC ENERGY" in line:
                electronic_energy = self.n_number_from_line(line, 0, 1)

            elif "CORE-CORE REPULSION" in line:
                core_repulsion = self.n_number_from_line(line, 0, 1)

        #total energy is the sum of electronic energy and core repulsion
        if electronic_energy is not None and core_repulsion is not None:
            self.energy = self.ev_to_au(electronic_energy + core_repulsion)
            self.log("NOTE: energies from semiempirical methods are not directly comparable to ab initio energies")

        else:
            self.log("Unable to find energy. Electronic energy: {0} Core-core repulsion: {1}".format(electronic_energy, core_repulsion))

    def extract_heat_of_formation(self, data, options={}):
        """Get heat of formation from log file and store it as
        self.heat_of_formation.

        :param data: log file contents
        :type data : str
        :param options: ignored
        :type options : dict
        """

        for line in data.split("\n"):
            if "HEAT OF FORMATION" in line:
                hof = self.n_number_from_line(line, 0, 1)
                self.heat_of_formation = self.kcalm_to_au(hof)

    def run(self, host="localhost", options={}):
        """Run MOPAC7 on the given host using the run_mopac7 script.

        :param host: name of host where job should execute
        :type host : str
        :param options: ignored
        :type options : dict
        """

        run_params = self.get_run_config(host)
        workdir = self.backend + "-" + str(uuid.uuid1()).replace('-', '')[:16]
        path = "{0}/{1}/".format(self.tmpdir, workdir)

        deck_hash = hashlib.sha1(self.deck).hexdigest()[:10]
        dat_file = "{0}.dat".format(deck_hash)
        abs_file = path + dat_file
        self.write_file(self.deck, abs_file, host)

        #N.B.: run_mopac7 does not like long paths!
        base_file = abs_file.split(".dat")[0]
        out_file = "{0}.out".format(base_file)
        rp = {"path" : path, "input" : base_file, "output" : out_file}
        cmd = run_params["cli"].format(**rp)

        stdout, returncode = self.execute(cmd, host, bash_shell=True)
        self.stdout = stdout
        self.logdata = self.read_file(out_file, host)

        errors = ["DUE TO PROGRAM BUG", "STOPPED TO AVOID WASTING TIME"]
        for e in errors:
            if e in self.logdata:
                self.runstate = "error"
                return
                
        self.extract_last_energy(self.logdata)
        self.extract_heat_of_formation(self.logdata)
        self.extract_geometry(self.logdata)

        self.runstate = "complete"

class Mopac7(cpinterface.MolecularCalculator):
    def __init__(self, *args, **kw):
        super(Mopac7, self).__init__(*args, **kw)
        
        self.methods = ["semiempirical:am1", "semiempirical:pm3",
                        "semiempirical:mindo/3", "semiempirical:mndo"]
        self.coordinate_choices = ["cartesian", "zmatrix"]
        self.references = ["RHF", "UHF"]

    def check_element_support(self, system, method):
        """Check that the chosen semiempirical method is parameterized for
        all the elements in the system. Supported elements are taken from
        section 3.5 of the Mopac7 manual. Unsupported elements will raise
        an exception.

        Note that MINDO/3 is supported only for certain *pairs* of elements,
        and this check may let bad pairs slip through because it is not
        pair-aware.

        :param system: molecular system
        :type system : geoprep.System
        :param method: name of method
        :type method : str
        :return: elements from system
        """
        
        emap = {"semiempirical:mndo" : ["H", "Li", "B", "C", "N", "O", "F",
                                        "Al", "Si", "P", "S", "Cl", "Zn",
                                        "Ge", "Br", "Sn", "I", "Hg", "Pb"],
                "semiempirical:am1" : ["H", "B", "C", "N", "O", "F", "Al",
                                       "Si", "P", "S", "Cl", "Zn", "Ge", "Br",
                                       "Sn", "I", "Hg"],
                "semiempirical:pm3" : ["H", "Be", "C", "N", "O", "F", "Mg",
                                       "Al", "Si", "P", "S", "Cl", "Zn", "Ga",
                                       "Ge", "As", "Se", "Br", "Cd", "In",
                                       "Sn", "Sb", "Te", "I", "Hg", "Tl", "Pb",
                                       "Bi"],
                "semiempirical:mindo/3" : ["H", "B", "C", "N", "O", "F", "Si",
                                           "P", "S", "Cl"]}
        
        elements = system.elements
        allowed = emap[method]
        for e in elements:
            if e not in allowed:
                raise ValueError("Element {0} not parameterized for {1}".format(repr(e), repr(method)))

        return elements

    def check_coordinates(self, coordinate_choice):
        if coordinate_choice in self.coordinate_choices:
            return coordinate_choice

        else:
            raise ValueError("Unrecognized coordinates option {0}".format(repr(coordinate_choice)))

    
    def create_geometry(self, system, options={}):
        """Create input geometry for a subsequent calculation.

        options:
         coordinates: "cartesian" or "zmatrix"

        :param system: molecular system data to convert to input geometry
        :type system : geoprep.System
        :param options: select coordinate system
        :type options : dict
        :return: a Mopac7 input with geometry specifications
        :rtype : str
        """

        defaults = {"coordinates" : "cartesian"}
        options = dict(defaults.items() + options.items())
        
        coord_choice = self.check_coordinates(options.get("coordinates"))

        if coord_choice == "cartesian":
            geometry = system.write("mopcrt")

        elif coord_choice == "zmatrix":
            geometry = system.write("mopin")

        return geometry

    def make_energy_job(self, system, method, options={}):
        """Create an input specification for a single point energy calculation.

        :param system: molecular system for energy calculation
        :type system : geoprep.System
        :param method: calculation method
        :type method : str
        :param options: additional keyword based control options
        :type options : dict
        :return: a Mopac7 single point energy calculation job
        :rtype : cpinterface.Job
        """
        
        system = self.fragment_to_system(system)

        if method.startswith("semiempirical"):
            return self.make_semiempirical_job(system, method, "ENERGY",
                                               options=options)
        else:
            raise ValueError("Mopac7 does not support {0}".format(method))

    def make_opt_job(self, system, method, options={}):
        """Create an input specification for a geometry optimization
        calculation. Optimization goal may be to find minimum geometry or
        to find a saddle point.

        options:
         goal: minimize or saddle

        :param system: molecular system for energy calculation
        :type system : geoprep.System
        :param method: calculation method
        :type method : str
        :param options: additional keyword based control options
        :type options : dict
        :return: a Mopac7 input for geometry optimization calculation
        :rtype : str
        """

        system = self.fragment_to_system(system)

        if method.startswith("semiempirical"):
            return self.make_semiempirical_job(system, method, "OPT",
                                               options=options)
        else:
            raise ValueError("Mopac7 does not support {0}".format(method))


    def make_semiempirical_job(self, system, method, runtyp,
                               options={}):
        """Create a semiempirical input specification for a calculation.
        Mopac7 supports MNDO, MINDO/3, AM1, and PM3 methods.

        See Chapter 2 of the Mopac7 manual for keyword details.

        :param system: molecular system for calculation
        :type system : geoprep.System
        :param method: a semiempirical calculation method
        :type method : str
        :param options: additional keyword based control options
        :type options : dict
        :return: a Mopac7 semiempirical job
        :rtype : Job
        """

        defaults = {"reference" : "rhf", "gnorm" : 0.0001, "precise" : True,
                    "let" : True, "scf_iterations" : 999, "geo_ok" : True,
                    "goal" : "minimize"}

        options = dict(defaults.items() + options.items())

        self.check_method(method)
        self.check_element_support(system, method)

        deck = self.create_geometry(system, options=options)
        semethod = method.split("semiempirical:")[-1].upper()        

        #MNDO is default method in Mopac7, so no keyword provided
        mmap = {"MNDO" : "" , "AM1" : "AM1",
                "MINDO/3" : "MINDO3", "PM3" : "PM3"}

        #RHF is default electronic reference in Mopac7, so no keyword provided
        rmap = {"RHF" : "", "UHF" : "UHF"}

        controls = []

        if system.spin > 6:
            raise ValueError("Spin of {0} too large: Mopac7 only works up to sextet".format(system.spin))

        else:
            spin_map = {1 : "SINGLET", 2 : "DOUBLET", 3 : "TRIPLET",
                        4 : "QUARTET", 5 : "QUINTET", 6 : "SEXTET"}
            spin_name = spin_map[system.spin]


        reference = rmap[options["reference"].upper()]
        self.check_electronic_reference(reference or options["reference"].upper())
        if system.spin > 1 and reference != "UHF":
            self.log("Forcing UHF for multiplicity {0}".format(system.spin))
            reference = "UHF"
            
        controls.append(reference)
        controls.append(spin_name)
        controls.append(mmap[semethod])
        controls.append("LARGE")

        #Default is to optimize geometry; "1SCF" gives a single-point energy
        if runtyp == "ENERGY":
            controls.append("1SCF")

        if options.get("gnorm"):
            controls.append("GNORM={0:.5f}".format(options.get("gnorm")))

        if options.get("precise"):
            controls.append("PRECISE")

        if options.get("let"):
            controls.append("LET")

        if options.get("geo_ok"):
            controls.append("GEO-OK")
        
        if options.get("scf_iterations"):
            controls.append("ITRY={0}".format(options.get("scf_iterations")))

        if system.charge != 0:
            controls.append("CHARGE={0}".format(system.charge))

        keywords = " ".join([c for c in controls if c])
        deck = deck.replace("PUT KEYWORDS HERE", keywords)

        job = Mopac7Job(deck=deck, system=system)

        return job

