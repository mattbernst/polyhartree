# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-

import hashlib
import os
import uuid
import shutil
import sys
import cpinterface

class PDynamoJob(cpinterface.Job):
    def __init__(self, *args, **kw):
        super(PDynamoJob, self).__init__(*args, **kw)
        self.backend = "pdynamo"

    def extract_last_energy(self, data, options={}):
        """Get last energy message from log file and store it as energy.

        @param data: log file contents
        @type data : str
        """

        electronic_energy = None
        nuclear_energy = None
        
        for line in data.split("\n"):
            if "Electronic Energy" in line and "Nuclear Energy" in line:
                electronic_energy = self.n_number_from_line(line, 0, 2)
                nuclear_energy = self.n_number_from_line(line, 1, 2)

        #total energy is the sum of electronic and nuclear (core) energy
        #units are already Hartree
        if electronic_energy is not None and nuclear_energy is not None:
            self.energy = electronic_energy + nuclear_energy
            self.log("NOTE: energies from semiempirical methods are not directly comparable to ab initio energies")

        else:
            self.log("Unable to find energy. Electronic energy: {0} Nuclear energy: {1}".format(electronic_energy, nuclear_energy))

    def extract_heat_of_formation(self, data, options={}):
        """Get heat of formation from log file and store it
        as heat_of_formation.

        @param data: log file contents
        @type data : str
        """

        for line in data.split("\n"):
            if "Heat of Formation:" in line:
                hof = self.n_number_from_line(line, 0, 1)
                self.heat_of_formation = self.kcalm_to_au(hof)
        
    def run_local(self, options={}):
        """The pdynamo-runner.py gets geometry from an .xyz file and everything
        else from command line arguments. This runner needs to copy the
        geometry and the runner to the working directory before executing the
        runner.
        """
        
        workdir = self.backend + "-" + str(uuid.uuid1()).replace('-', '')[:16]
        path = "/tmp/{0}/".format(workdir)
        os.makedirs(path)
        
        #write .xyz geometry file to working directory
        xyzdata = self.system.write("xyz")
        xyzfile = "{0}.xyz".format(workdir)
        xyzfull = "{0}/{1}".format(path, xyzfile)
        with open(xyzfull, "w") as geofile:
            geofile.write(xyzdata)

        #add geometry file specification to command line
        log_file = path + xyzfile.replace(".xyz", ".log")
        self.deck += " --xyzfile={0} --outfile={1}".format(xyzfull, log_file)

        #copy runner to working directory
        try:
            here = os.environ["POLYHARTREE_ROOT"]
            src = "{0}/pdynamo-runner.py".format(here)
            dst = path + "pdynamo-runner.py"
            shutil.copy(src, dst)
        except KeyError:
            msg = "You need to set the environment variable POLYHARTREE_ROOT to point to the top level polyhartree directory, e.g. export POLYHARTREE_ROOT=/home/steve/polyhartree\n"
            sys.stderr.write(msg)
            sys.exit(1)
            
        cmd = self.deck
        
        stdout, returncode = self.execute(cmd)
        
        self.stdout = stdout

        with open(log_file, "r") as lfile:
            self.logdata = lfile.read()
            self.extract_last_energy(self.logdata)
            self.extract_heat_of_formation(self.logdata)

        self.runstate = "complete"

class PDynamo(cpinterface.MolecularCalculator):
    def __init__(self, *args, **kw):
        super(PDynamo, self).__init__(*args, **kw)
        
        self.methods = ["semiempirical:am1", "semiempirical:pm3",
                        "semiempirical:pm6", "semiempirical:mndo",
                        "semiempirical:rm1", "semiempirical:pddg/pm3",
                        "semiempirical:pddg/mndo", "semiempirical:am1/d-phot"]
        self.coordinate_choices = ["cartesian"]
        self.references = ["RHF", "UHF"]

    def check_coordinates(self, coordinate_choice):
        if coordinate_choice in self.coordinate_choices:
            return coordinate_choice

        else:
            raise ValueError("Unrecognized coordinates option {0}".format(repr(coordinate_choice)))

    
    def create_geometry(self, system, options={}):
        """This method doesn't really apply to pDynamo calculations. The job
        runner needs to create an .xyz geometry file since there is no proper
        input deck for pdynamo calculations, only an .xyz file plus a large
        number of command line arguments.
        """

        return ""

    def make_energy_job(self, system, method, options={}):
        """Create an input specification for a single point energy calculation.

        @param system: molecular system for energy calculation
        @type system : cinfony molecule
        @param method: calculation method
        @type method : str
        @return: a pDynamo single point energy calculation job
        @rtype : cpinterface.Job
        """

        if method.startswith("semiempirical"):
            return self.make_semiempirical_job(system, method, "ENERGY",
                                               options=options)
        else:
            raise ValueError("pDynamo 7 does not support {0}".format(method))

    def make_semiempirical_job(self, system, method, runtyp,
                               options={}):
        """Create a semiempirical input specification for a calculation.
        pDynamo supports AM1, MNDO, PM3, PM6, RM1, AM1/d-PhoT, PDDG/PM3,
        PDDG/MNDO

        @param system: molecular system for energy calculation
        @type system : cinfony molecule
        @param method: a semiempirical calculation method
        @type method : str
        @return: a PDynamo semiempirical job
        @rtype : Job
        """

        defaults = {"reference" : "rhf"}
        options = dict(defaults.items() + options.items())

        self.check_method(method)

        semethod = method.split("semiempirical:")[-1].upper()        

        mmap = {"MNDO" : "mndo" , "AM1" : "am1", "PM3" : "pm3", "PM6" : "pm6",
                "AM1/D-PHOT" : "am1dphot", "PDDG/PM3" : "pddgpm3",
                "PDDG/MNDO" : "pddgmndo", "RM1" : "rm1"}

        rmap = {"RHF" : "RHF", "UHF" : "UHF"}

        controls = []

        reference = rmap[options["reference"].upper()]
        self.check_electronic_reference(reference or options["reference"].upper())
        if system.spin > 1 and reference == "RHF":
            self.log("Forcing UHF for multiplicity {0}".format(system.spin))
            reference = "UHF"
            
        controls.append("--restricted={0}".format(reference))
        controls.append("--multiplicity={0}".format(system.spin))
        mname = "semiempirical:{0}".format(mmap[semethod])
        controls.append("--method={0}".format(mname))

        if runtyp == "ENERGY":
            controls.append("--jobtype=energy")

        controls.append("--charge={0}".format(system.charge))

        #The "deck" here is just a series of command line arguments for
        #pdynamo-runner.py. In the actual job runner we will add the name
        #of the geometry file, create a shell script, and copy the runner.
        args = " ".join([c for c in controls if c])
        deck = "python pdynamo-runner.py " + args

        job = PDynamoJob(deck=deck, system=system)

        return job

