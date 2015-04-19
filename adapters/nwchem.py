# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-

import hashlib
import os
import uuid
import cpinterface

class NWChemJob(cpinterface.Job):
    def __init__(self, *args, **kw):
        super(NWChemJob, self).__init__(*args, **kw)
        self.backend = "nwchem"

    def extract_last_energy(self, data, options={}):
        """Get last energy message from log file and store it as energy.

        @param data: log file contents
        @type data : str
        @param options: ignored
        @type options : dict
        """

        for line in data.split("\n"):
            if "Total SCF energy" in line:
                energy = self.n_number_from_line(line, 0, 1)

                #units are already Hartree
                self.energy = energy

    def run_local(self, options={}):
        """Run a NWChem job using mpirun, on the local host.

        @param options: ignored
        @type options : dict
        """
        
        workdir = self.backend + "-" + str(uuid.uuid1()).replace('-', '')[:16]
        path = "/tmp/{0}/".format(workdir)
        os.makedirs(path)

        deck_hash = hashlib.sha1(self.deck).hexdigest()[:10]
        dat_file = "{0}.nw".format(deck_hash)
        abs_file = path + dat_file
        log_file = abs_file.replace(".nw", ".log")
        with open(abs_file, "w") as deckfile:
            deckfile.write(self.deck)

        cmd = "cd {0} && mpirun -np 1 nwchem {1} &> {2}".format(path, dat_file,
                                                                log_file)
        
        stdout, returncode = self.execute(cmd, cwd=path, bash_shell=True)
        self.stdout = stdout

        with open(log_file, "r") as lfile:
            self.logdata = lfile.read()

            if "There is an error in the input" in self.logdata:
                self.runstate = "error"
                
            else:
                self.extract_last_energy(self.logdata)
                self.runstate = "complete"

class NWChem(cpinterface.MolecularCalculator):
    def __init__(self, *args, **kw):
        super(NWChem, self).__init__(*args, **kw)
        
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
        @return: a NWChem input with geometry specifications
        @rtype : str
        """

        defaults = {"coordinates" : "cartesian"}
        options = dict(defaults.items() + options.items())
        symmetry = options.get("symmetry")
        
        coord_choice = self.check_coordinates(options.get("coordinates"))

        if coord_choice == "cartesian":
            coordinates = self.make_tagged_cartesian_geometry(system,
                                                              options=options)
            head = "geometry units angstroms print xyz autosym"
            directives = [head]
            if symmetry:
                directives.append("symmetry {0}".format(symmetry))
            geometry = self.make_control_block(directives + coordinates)

        elif coord_choice == "zmatrix":
            raise ValueError("NWChem zmatrix input currently unsupported")

        return geometry

    def make_control_block(self, controls, indent="  "):
        """Make a textual control block for NWChem. For example, to modify
        SCF parameters, pass in ["scf", "thresh 1e-6" "maxiter 200"] as
        controls and a block like this will be generated:

        scf
          thresh 1e-6
          maxiter 200
        end

        @param controls: block name followed by parameters to insert in block
        @type controls : list
        @param indent: whitespace to prepend to block parameters
        @type indent : str
        @return: formatted control block
        @rtype : str
        """

        lines = [controls[0]]
        for control in controls[1:]:
            line = indent + control
            lines.append(line)

        lines.append("end")
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
        @return: a NWChem input for single point energy calculation
        @rtype : str
        """

        system = self.fragment_to_system(system)

        self.check_method(method)

        if method.startswith("hf"):
            return self.make_hf_job(system, method, "ENERGY", options=options)

        else:
            raise ValueError("NWChem does not currently support {0}".method)

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
        bsd = self.get_basis_data(system, {"basis_format" : "nwchem"})
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
                el_size = len(element)
                basis_text = data[basis_name][element]

                #add basis data with comment e.g. "H 3-21G"
                #each element symbol gets replaced by tag, e.g. Cl -> Cl0
                comment = "#{0} {1}".format(element, basis_name)
                basis_lines = [comment]
                for line in basis_text.split("\n")[1:]:
                    if line[:el_size] == element:
                        suffix = line[el_size:]
                        entry = basis_alias + suffix
                    else:
                        entry = line
                        #remove premature END terminators
                        if entry == "END":
                            continue
                    basis_lines.append(entry)
                    
                block = "\n".join(basis_lines)
                
                basis_blocks.append(block)

        blocks = "\n".join(basis_blocks)
        head = "basis {0}".format(form)
        formatted = self.make_control_block([head] + blocks.split("\n"))

        #set basis tag for each atom
        system.set_properties(property_name, range(len(system.atoms)),
                              basis_names)

        r = {"basis_data" : formatted}
        
        return r

    def make_hf_job(self, system, method, runtyp, options={}):
        """Create an input specification for a Hartree-Fock calculation.
        
        NWChem supports RHF restricted Hartree-Fock,
        UHF unrestricted Hartree-Fock, and
        ROHF restricted open shell Hartree-Fock

        @param system: molecular system for calculation
        @type system : geoprep.System
        @param method: a HF calculation method
        @param options: additional keyword based control options
        @type options : dict
        @type method : str
        @return: a NWChem HF job
        @rtype : Job
        """

        defaults = {"reference" : "rhf", "scf_iterations" : 999,
                    "basis_tag_name" : "basis_tag",
                    "property_name" : "basis_tag"}
        options = dict(defaults.items() + options.items())

        self.check_method(method)

        bd = self.prepare_basis_data(system, options=options)
        geometry = self.create_geometry(system, options=options)
        
        method = method.split("hf:")[-1].upper()

        reference = options.get("reference").upper()
        self.check_electronic_reference(reference or options["reference"].upper())
        if system.spin > 1 and reference == "RHF":
            self.log("Forcing UHF for multiplicity {0}".format(system.spin))
            reference = "UHF"

        scf = ["scf", reference, "nopen {0}".format(system.spin - 1),
               "maxiter {0}".format(options.get("scf_iterations"))]
        scf = self.make_control_block(scf)
        charge = "charge {0}".format(system.charge)

        deck = ["start molecule",
                """title "{0}" """.format(system.title),
                geometry,
                scf,
                charge,
                bd["basis_data"],
                "task scf energy"]
        
        deck = "\n\n".join(deck)
        
        job = NWChemJob(deck=deck, system=system)
        return job
