# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-

import hashlib
import uuid

import cpinterface
import sharedutilities

class Psi4Job(cpinterface.Job):
    def __init__(self, *args, **kw):
        super(Psi4Job, self).__init__(*args, **kw)
        self.backend = "psi4"

    def extract_last_energy(self, data, options={}):
        """Get last energy message from log file and store it as self.energy.

        @param data: log file contents
        @type data : str
        @param options: ignored
        @type options : dict
        """

        for line in data.split("\n"):
            if "Final Energy" in line:
                energy = self.n_number_from_line(line, 0, 1)

                #units are already Hartree
                self.energy = energy

    def extract_geometry(self, data, options={}):
        """Get last geometry found in log file and store it as self.geometry.

        @param data: log file contents
        @type data : str
        @param options: ignored
        @type options : dict
        """

        u = sharedutilities.Utility()
        initial_geometry = []
        geometries = []
        init = False
        # initial geometry starts after XYZ format geometry
        for line in data.split("\n"):
            if "Geometry (in Angstrom)" in line:
                init = True

            elif init:
                #run block follows geometry
                if "Running in" in line:
                    break

                coords = u.numericize(line)
                if sum([1 for c in coords if type(c) == float]) == 3:
                    #got a coordinate line: symbol, x, y, z
                    initial_geometry.append(coords)

        if geometries:
            self.geometry = geometries[-1]
        else:
            self.geometry = initial_geometry



    def run(self, host="localhost", options={}):
        """Run a Psi4 job using psi script, on the local host.

        @param host: name of host where job should execute
        @type host : str
        @param options: ignored
        @type options : dict
        """

        run_params = self.get_run_config(host)
        workdir = self.backend + "-" + str(uuid.uuid1()).replace('-', '')[:16]
        path = "{0}/{1}/".format(self.tmpdir, workdir)

        deck_hash = hashlib.sha1(self.deck).hexdigest()[:10]
        dat_file = "{0}.dat".format(deck_hash)
        abs_file = path + dat_file
        log_file = abs_file.replace(".dat", ".log")
        self.write_file(self.deck, abs_file, host)

        rp = {"path" : path, "input" : dat_file, "output" : log_file,
              "ncores" : run_params["cores"]}
        cmd = run_params["cli"].format(**rp)
                
        stdout, returncode = self.execute(cmd, host, cwd=path, bash_shell=True)
        self.stdout = stdout

        self.logdata = self.read_file(log_file, host)

        errors = ["PsiException:", "Error:"]

        for e in errors:
            if e in self.logdata:
                self.runstate = "error"

        if self.runstate != "error":
            self.extract_last_energy(self.logdata)
            self.extract_geometry(self.logdata)
            self.runstate = "complete"

class Psi4(cpinterface.MolecularCalculator):
    def __init__(self, *args, **kw):
        super(Psi4, self).__init__(*args, **kw)
        
        self.methods = ["hf:rhf", "hf:uhf", "hf:rohf"]
        self.coordinate_choices = ["cartesian", "zmatrix"]
        self.references = ["RHF", "ROHF", "UHF"]

    def make_tagged_cartesian_geometry(self, system, options={}):
        """Create a geometry block using tagged atom names.

        options:
         property_name: name of tag group to use from atom_properties

        @param system: molecular system data to convert to input geometry
        @type system : geoprep.System
        @return: atom coordinates named by atom tags
        @rtype : list
        """

        tag_name = options["property_name"]
        tags = system.atom_properties(tag_name)

        entries = []
        for j, tag in enumerate(tags):
            atom = system.atoms[j]
            entry = "{}\t{:.6f}\t{:.6f}\t{:.6f}".format(tag, atom.coords[0],
                                                        atom.coords[1],
                                                        atom.coords[2])
            entries.append(entry)

        return entries

    def create_geometry(self, system, options={}):
        """Create input geometry for a subsequent calculation.

        options:
         coordinates: "cartesian" or "zmatrix"
         symmetry: optional symmetry group

        @param system: molecular system data to convert to input geometry
        @type system : geoprep.System
        @param options: select coordinate system
        @type options : dict
        @return: a Psi4 input with geometry specifications
        @rtype : str
        """

        defaults = {"coordinates" : "cartesian"}
        options = dict(defaults.items() + options.items())
        symmetry = options.get("symmetry")
        
        coord_choice = self.check_coordinates(options.get("coordinates"))

        if coord_choice == "cartesian":
            coordinates = self.make_tagged_cartesian_geometry(system,
                                                              options=options)
            head = "molecule"
            #psi4 reads charge and multiplicity from geometry block
            chargemult = "{0} {1}".format(system.charge, system.spin)
            directives = [head, "#charge multiplicity", chargemult]
            if symmetry:
                directives.append("symmetry {0}".format(symmetry))
            geometry = self.make_control_block(directives + coordinates)

        elif coord_choice == "zmatrix":
            raise ValueError("Psi4 zmatrix input currently unsupported")

        return geometry

    def make_control_block(self, controls, indent="  "):
        """Make a textual control block for Psi4. For example, to modify
        assign basis parameters, pass in
        ["basis", "assign H1 sto-3g" "assign C1 sto-3g"] as
        controls and a block like this will be generated:

        basis {
          assign H1 sto-3g
          assign C1 sto-3G
        }

        @param controls: block name followed by parameters to insert in block
        @type controls : list
        @param indent: whitespace to prepend to block parameters
        @type indent : str
        @return: formatted control block
        @rtype : str
        """

        lines = [controls[0] + " {"]
        for control in controls[1:]:
            line = indent + control
            lines.append(line)

        lines.append("}")
        block = "\n".join(lines)
        return block

    def make_energy_job(self, system, method, options={}):
        """Create an input specification for a single point energy calculation.

        @param system: molecular system for energy calculation
        @type system : geoprep.System
        @param method: calculation method
        @type method : str
        @param options: additional keyword based control options
        @type options : dict
        @return: a Psi4 input for single point energy calculation
        @rtype : str
        """

        system = self.fragment_to_system(system)

        self.check_method(method)

        if method.startswith("hf"):
            return self.make_hf_job(system, method, "ENERGY", options=options)

        else:
            raise ValueError("Psi4 does not currently support {0}".method)

    def prepare_basis_data(self, system, options={}):
        """Prepare basis set data for use: basis assignment, inline
        basis specifications, and information for spherical/cartesian flag.

        Also assign basis tags to atoms as atom properties.

        @param system: molecular system for parameterization
        @type system : geoprep.System
        @param options: additional keyword based control options
        @type options : dict
        @return: basis set control block
        @rtype : dict
        """

        property_name = options["basis_tag_name"]
        ubb = {}
        basis_names = []
        basis_blocks = []
        assignments = []
        bsd = self.get_basis_data(system, {"basis_format" : "g94"})
        data = bsd["data"]

        form = bsd["spherical_or_cartesian"]

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
            
            if add_new_block:
                assignment = "assign {0} {1}".format(basis_alias, basis_alias)
                assignments.append(assignment)
                basis_text = data[basis_name][element]

                #add basis data with comment e.g. "H 3-21G"
                #also indicate cartesian/spherical form
                comment = "#{0} {1}".format(element, basis_name)
                designation = "[ {0} ]".format(basis_alias)
                basis_lines = [comment, designation, form, basis_text]
                                    
                block = "\n".join(basis_lines)
                
                basis_blocks.append(block)

        blocks = "\n".join(basis_blocks)
        formatted = blocks

        data = self.make_control_block(["basis"] + assignments +
                                       blocks.split("\n"))

        #set basis tag for each atom
        system.set_properties(property_name, range(len(system.atoms)),
                              basis_names)

        r = {"basis_data" : data}
        
        return r

    def make_hf_job(self, system, method, runtyp, options={}):
        """Create an input specification for a Hartree-Fock calculation.
        
        Psi4 supports RHF restricted Hartree-Fock,
        UHF unrestricted Hartree-Fock, and
        ROHF restricted open shell Hartree-Fock

        options:
         scf_type: How the electron repulsion integrals are calculated. The
          psi4 default is df, density fitting, but that is not directly
          comparable to methods used by GAMESS or NWChem. Instead polyhartree
          uses pk as default. Possible options are pk, out_of_core, direct,
          df, and cd.

        @param system: molecular system for calculation
        @type system : geoprep.System
        @param method: a HF calculation method
        @param options: additional keyword based control options
        @type options : dict
        @type method : str
        @return: a Psi4 HF job
        @rtype : Job
        """

        defaults = {"scf_iterations" : 999,
                    "basis_tag_name" : "basis_tag",
                    "property_name" : "basis_tag",
                    "scf_type" : "pk"}
        options = dict(defaults.items() + options.items())

        self.check_method(method)

        bd = self.prepare_basis_data(system, options=options)
        geometry = self.create_geometry(system, options=options)
        
        reference = method.split("hf:")[-1].upper()

        self.check_electronic_reference(reference or options["reference"].upper())
        if system.spin > 1 and reference == "RHF":
            self.log("Forcing UHF for multiplicity {0}".format(system.spin))
            reference = "UHF"

        gb = ["set globals",
              "scf_type {0}".format(options.get("scf_type")),
              "reference {0}".format(reference),
              "maxiter {0}".format(options.get("scf_iterations"))]
        gb = self.make_control_block(gb)

        deck = ["# {0}".format(system.title),
                geometry,
                gb,
                bd["basis_data"],
                "energy('scf')"]
        
        deck = "\n\n".join(deck)
        
        job = Psi4Job(deck=deck, system=system)
        return job
