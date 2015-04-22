# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-

import hashlib
import uuid
import cpinterface

class GAMESSUSJob(cpinterface.Job):
    def __init__(self, *args, **kw):
        super(GAMESSUSJob, self).__init__(*args, **kw)
        self.backend = "gamess_us"

    def extract_last_energy(self, data, options={}):
        """Get last energy message from log file and store it as energy.

        @param data: log file contents
        @type data : str
        @param options: ignored
        @type options : dict
        """

        for line in data.split("\n"):
            if "FINAL " in line and "ENERGY IS" in line:
                energy = self.n_number_from_line(line, 0, 2)

                #units are already Hartree
                self.energy = energy

                #Note: slight differences from other Mopac implementations
                #because GAMESS-US, as of release 1 MAY 2013 (R1),
                #treats the Hartree as 27.211652 eV whereas current accepted
                #value is 27.211385 eV (see value TOHART in GAMESS mpcint.src)
                if self.extras.get("semiempirical"):
                    self.log_once("NOTE: energies from semiempirical methods are not directly comparable to ab initio energies")

    def extract_heat_of_formation(self, data, options={}):
        """Get heat of formation from log file and store it
        as heat_of_formation.

        @param data: log file contents
        @type data : str
        @param options: ignored
        @type options : dict
        """

        for line in data.split("\n"):
            if "HEAT OF FORMATION IS" in line:
                hof = self.n_number_from_line(line, 0, 1)
                self.heat_of_formation = self.kcalm_to_au(hof)

    def run(self, host="localhost", options={}):
        """Run a GAMESS-US job on the given host.

        @param host: name of host where job should execute
        @type host : str
        @param options: ignored
        @type options : dict
        """

        if host != "localhost":
            raise NotImplementedError("Remote job execution not yet ready")

        run_params = self.get_run_config(host)
        workdir = self.backend + "-" + str(uuid.uuid1()).replace('-', '')[:16]
        path = "/tmp/{0}/".format(workdir)

        deck_hash = hashlib.sha1(self.deck).hexdigest()[:10]
        dat_file = "{0}.inp".format(deck_hash)
        abs_file = path + dat_file
        log_file = abs_file.replace(".inp", ".log")
        self.write_file(self.deck, abs_file, host)

        rp = {"path" : path, "input" : deck_hash, "output" : log_file,
              "ncores" : run_params["cores"]}
        cmd = run_params["cli"].format(**rp)
        
        stdout, returncode = self.execute(cmd, cwd=path, bash_shell=True)
        self.stdout = stdout

        errors = ["FATAL ERROR", "TYPING ERROR",
                  "GAMESS TERMINATED -ABNORMALLY-"]
        self.logdata = self.read_file(log_file, host)
 
        logupper = self.logdata.upper()
        for e in errors:
            if e in logupper:
                self.runstate = "error"

        if self.runstate != "error":
            self.extract_last_energy(self.logdata)
            self.extract_heat_of_formation(self.logdata)
            self.runstate = "complete"

class GAMESSUS(cpinterface.MolecularCalculator):
    def __init__(self, *args, **kw):
        super(GAMESSUS, self).__init__(*args, **kw)
        
        self.methods = ["semiempirical:am1", "semiempirical:pm3",
                        "semiempirical:rm1", "semiempirical:mndo",
                        "hf:rhf", "hf:uhf", "hf:rohf"]
        self.coordinate_choices = ["cartesian", "zmatrix"]
        self.references = ["RHF", "ROHF", "UHF"]

    def check_element_support(self, system, method):
        """Check that the chosen semiempirical method is parameterized for
        all the elements in the system. Supported elements are taken from
        the GAMESS manual, "MOPAC Calculations Within GAMESS",
        http://www.msg.ameslab.gov/gamess/GAMESS_Manual/refs.pdf

        @param system: molecular system
        @type system : geoprep.System
        @param method: name of method
        @type method : str
        @return: elements from system
        """

        emap = {"semiempirical:pm3" : ["H", "Li", "Be", "C", "N", "O", "F",
                                       "Na", "Mg", "Al", "Si", "P", "S", "Cl",
                                       "K", "Ca", "Zn", "Ga", "Ge", "As", "Se",
                                       "Br", "Cd", "In", "Sn", "Sb", "Te",
                                       "I", "Hg", "Tl", "Pb", "Bi"],
                "semiempirical:am1" : ["H", "B", "C", "N", "O", "F", "Na",
                                       "Mg", "Al", "Si", "P", "S", "Cl",
                                       "K", "Ca", "Zn", "Ge", "Br", "Sn", "I",
                                       "Hg"],
                "semiempirical:mndo" : ["H", "Li", "B", "C", "N", "O", "F",
                                        "Al", "Si", "P", "S", "Cl", "Zn",
                                        "Ge", "Br", "Sn", "I", "Hg", "Pb"]}

        #RM1 has its own parameters for C, N, O, F, P, S, Cl, Br, I;
        #otherwise AM1 parameters are used
        emap["semiempirical:rm1"] = emap["semiempirical:am1"]

        elements = system.elements
        allowed = emap[method]
        for e in elements:
            if e not in allowed:
                raise ValueError("Element {0} not parameterized for {1}".format(repr(e), repr(method)))

        return elements
    
    def create_geometry(self, system, options={}):
        """Create input geometry for a subsequent calculation.

        options:
         coordinates: "cartesian" or "zmatrix"
         symmetry: optional symmetry group

        @param system: molecular system data to convert to input geometry
        @type system : geoprep.System
        @param options: select coordinate system
        @type options : dict
        @return: a GAMESS-US input with geometry specifications
        @rtype : str
        """

        defaults = {"coordinates" : "cartesian"}
        options = dict(defaults.items() + options.items())
        symmetry = options.get("symmetry")
        
        coord_choice = self.check_coordinates(options.get("coordinates"))

        if coord_choice == "cartesian":
            geometry = system.write("inp")

            #change now-obsolete CART coordinate designation to
            #synonymous UNIQUE
            geometry=geometry.replace("COORD=CART", "COORD=UNIQUE")

            #set explicit symmetry
            if symmetry:
                lines = geometry.split("\n")
                lines[4] = symmetry
                geometry = "\n".join(lines)

        elif coord_choice == "zmatrix":
            raise ValueError("GAMESS-US zmatrix input currently unsupported")

        #add title
        geometry = geometry.replace(" $DATA\n\n",
                                    " $DATA\n" + system.title + "\n")

        return geometry

    def make_energy_job(self, system, method, options={}):
        """Create an input specification for a single point energy calculation.

        @param system: molecular system for energy calculation
        @type system : geoprep.System
        @param method: calculation method
        @type method : str
        @param options: additional keyword based control options
        @type options : dict
        @return: a GAMESS-US input for single point energy calculation
        @rtype : str
        """

        system = self.fragment_to_system(system)

        self.check_method(method)
        if method.startswith("semiempirical"):
            try:
                options["extras"]["semiempirical"] = True
            except KeyError:
                options["extras"] = {"semiempirical" : True}
                
            return self.make_semiempirical_job(system, method, "ENERGY",
                                               options=options)

        elif method.startswith("hf"):
            return self.make_hf_job(system, method, "ENERGY", options=options)

        else:
            raise ValueError("GAMESS-US does not currently support {0}".method)

    def reformat_long_line(self, data, start_marker, end_marker, maxlen=70):
        """GAMESS-US does not read lines longer than 80 characters. This
        method breaks up long directives into multiple lines.

        @param data: text of input deck that needs fixing
        @type data : str
        @param start_marker: token that marks fixup start, e.g. "$CONTRL"
        @type start_marker : str
        @param end_marker: token that marks fixup end, e.g. "$END"
        @type end_marker : str
        @param maxlen: maximum line length to generate
        @type maxlen : int
        @return: (shortened_line_list, region_start_index, region_end_index)
        @rtype : tuple
        """
            
        begin_index = data.find(start_marker)
        end_index = data.find(end_marker) + len(end_marker)
        inner = data[begin_index : end_index]

        pieces = []
        current_split = 0
        last_split = 0
        while inner:
            #stop when remainder fits wholly
            if len(inner) < maxlen:
                pieces.append(inner)
                break

            #look for largest split that isn't too large
            current_split = inner.find(" ", last_split + 1)
            if current_split >= maxlen or current_split == -1:
                t = inner[:last_split]
                pieces.append(t)
                inner = inner[last_split:]
                current_split = 0
                last_split = 0

            else:
                last_split = current_split

        return pieces, begin_index, end_index

    def make_semiempirical_job(self, system, method, runtyp, options={}):
        """Create a semiempirical input specification for a calculation.
        GAMESS-US supports MNDO, RM1, AM1, and PM3 methods.

        @param system: molecular system for calculation
        @type system : geoprep.System
        @param method: a semiempirical calculation method
        @type method : str
        @param options: additional keyword based control options
        @type options : dict
        @return: a GAMESS-US semiempirical job
        @rtype : Job
        """

        defaults = {"reference" : "rhf", "scf_iterations" : 200}
        options = dict(defaults.items() + options.items())

        self.check_method(method)
        self.check_element_support(system, method)

        deck = self.create_geometry(system, options=options)
        semethod = method.split("semiempirical:")[-1].upper()

        reference = options.get("reference").upper()
        self.check_electronic_reference(reference or options["reference"].upper())
        if system.spin > 1 and reference == "RHF":
            self.log("Forcing UHF for multiplicity {0}".format(system.spin))
            reference = "UHF"

        control_options = {"runtyp" : runtyp, "mult" : system.spin,
                           "icharg" : system.charge, "reference" : reference,
                           "scf_iterations" : options.get("scf_iterations")}
            
        contrl = "SCFTYP={reference} RUNTYP={runtyp} ICHARG={icharg} MULT={mult} MAXIT={scf_iterations} ".format(**control_options)
        deck = deck.replace("$CONTRL ", "$CONTRL " + contrl)

        pieces, begin, end = self.reformat_long_line(deck, "$CONTRL", "$END")
        new_contrl = "\n ".join(pieces)

        deck = deck[:begin] + new_contrl + deck[end:]
        basis = " $BASIS GBASIS={0} $END".format(semethod)
        lines = deck.split("\n")
        joblines = []
        inserted = False
        for line in lines:
            joblines.append(line)
            if "$END" in line and not inserted:
                inserted = True
                joblines.append(basis)

        deck = "\n".join(joblines)

        job = GAMESSUSJob(deck=deck, system=system,
                          extras=options.get("extras", {}))
        return job

    def prepare_basis_data(self, system, options={}):
        """Prepare basis set data for use: basis assignment array, inline
        basis specifications, and information for ISPHER flag.

        @param system: molecular system for parameterization
        @type system : geoprep.System
        @param options: additional keyword based control options
        @type options : dict
        @return: ISPHER flag data, assignment array, basis data, comments
        @rtype : dict
        """

        property_name = options["basis_tag_name"]
        ubb = {}
        basis_names = []
        basis_blocks = []
        comments = []
        bsd = self.get_basis_data(system, {"basis_format" : "gamess-us"})
        data = bsd["data"]
        if bsd["spherical_or_cartesian"] == "spherical":
            ispher = 1
        else:
            ispher = 0

        symbols = system.atom_properties("symbols")
        for j, basis_name in enumerate(system.atom_properties("basis_name")):
            add_new_block = False
            element = symbols[j]
            if element not in ubb:
                ubb[element] = []
                
            basis_element = "{0}:{1}".format(basis_name, element)

            #this (basis, element) pair is new, so add it
            if basis_element not in ubb[element]:
                ubb[element].append(basis_element)
                add_new_block = True

            basis_index = ubb[element].index(basis_element)
            basis_alias = "{0}{1}".format(element, basis_index)
            basis_names.append(basis_alias)
            comment = "{0} {1}".format(element, basis_name)
            comments.append(comment)
            
            if add_new_block:
                basis_text = data[basis_name][element]
                basis_text = "\n".join(basis_text.split("\n")[1:])
                
                tpl = "!\t{comment}\n ${name}\n{data}\n\n $END"
                block = tpl.format(**{"comment" : comment,
                                      "name" : basis_alias,
                                      "data" : basis_text})
                basis_blocks.append(block)

        blocks = "\n".join(basis_blocks)

        indexes = ", ".join(basis_names)
        bn = " $BASIS basnam(1)={0} $END".format(indexes)
        pieces, begin, end = self.reformat_long_line(bn, " $BASIS", "$END")
        new_bn = "\n ".join(pieces)

        r = {"ispher" : ispher, "basis_control" : new_bn,
             "basis_data" : blocks, "comments" : comments}

        #set basis tag for each atom
        system.set_properties(property_name, range(len(system.atoms)),
                              basis_names)

        return r

    def make_hf_job(self, system, method, runtyp, options={}):
        """Create an input specification for a Hartree-Fock calculation.
        
        GAMESS-US supports RHF restricted Hartree-Fock,
        UHF unrestricted Hartree-Fock, and
        ROHF restricted open shell Hartree-Fock

        @param system: molecular system for calculation
        @type system : geoprep.System
        @param method: a HF calculation method
        @param options: additional keyword based control options
        @type options : dict
        @type method : str
        @return: a GAMESS-US HF job
        @rtype : Job
        """

        defaults = {"reference" : "rhf", "scf_iterations" : 200,
                    "basis_tag_name" : "basis_tag"}
        options = dict(defaults.items() + options.items())

        self.check_method(method)

        deck = self.create_geometry(system, options=options)
        method = method.split("hf:")[-1].upper()

        reference = options.get("reference").upper()
        self.check_electronic_reference(reference or options["reference"].upper())
        if system.spin > 1 and reference == "RHF":
            self.log("Forcing UHF for multiplicity {0}".format(system.spin))
            reference = "UHF"

        bd = self.prepare_basis_data(system, options=options)

        control_options = {"runtyp" : runtyp, "mult" : system.spin,
                           "icharg" : system.charge, "reference" : reference,
                           "scf_iterations" : options.get("scf_iterations"),
                           "ispher" : bd["ispher"]}
            
        contrl = "SCFTYP={reference} RUNTYP={runtyp} ICHARG={icharg} ISPHER={ispher} MULT={mult} MAXIT={scf_iterations} ".format(**control_options)
        deck = deck.replace("$CONTRL ", "$CONTRL " + contrl)

        pieces, begin, end = self.reformat_long_line(deck, " $CONTRL", "$END")
        new_contrl = "\n ".join(pieces)

        deck = deck[:begin] + new_contrl + deck[end:]
        
        basis = bd["basis_control"]
        lines = deck.split("\n")
        joblines = []
        inserted = False
        for line in lines:
            joblines.append(line)
            if "$END" in line and not inserted:
                inserted = True
                joblines.append(basis)

        deck = "\n".join(joblines)
        deck += "\n!Basis set assignments:\n"
        deck += "\n".join(["!\t" + x for x in bd["comments"]]) + "\n\n"
        deck += bd["basis_data"]

        job = GAMESSUSJob(deck=deck, system=system)
        return job
