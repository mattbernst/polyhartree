# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import os
import uuid

import geoprep
import cpinterface
import sharedutilities

class PDynamoJob(cpinterface.Job):
    def __init__(self, *args, **kw):
        super(PDynamoJob, self).__init__(*args, **kw)
        self.backend = "pdynamo"

    def extract_last_energy(self, data, options={}):
        """Get last energy message from log file and store it as self.energy.

        :param data: log file contents
        :type data : str
        :param options: ignored
        :type options : dict
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
        as self.heat_of_formation.

        :param data: log file contents
        :type data : str
        :param options: ignored
        :type options : dict
        """

        scanning = False
        for line in data.split("\n"):
            if "Heat of Formation" in line:
                scanning = True
                
            elif scanning:
                try:
                    hof = self.n_number_from_line(line, 0, 1)
                except ValueError:
                    continue
                    
                self.heat_of_formation = self.kcalm_to_au(hof)
                break

    def line_to_geometry(self, line):
        """Extract lines fitting pattern element_symbol float float float
        and interpret them as [element, x, y, z].

        e.g. WILL match
         O   0.0133992460  -0.0168294554   0.0087021955

        WILL NOT match
         1  O            8.0    -0.0002646     0.0003329    -0.0001721

        :param line: a line of data from log file
        :type line : str
        :return: [element, x, y, z] or []
        :rtype : list
        """

        entry = []
        u = sharedutilities.Utility()
        n = u.numericize(line)
        pattern = [str, float, float, float]
        if [type(k) for k in n] == pattern:
            if n[0] in sharedutilities.ELEMENTS:
                entry = [n[0], n[1], n[2], n[3]]

        return entry

    def extract_geometry(self, data, options={}):
        """Get last geometry found in log file and store it as self.geometry.
        If there are multiple geometries from e.g. an optimization run, they
        will go into self.geometry_history.

        :param data: log file contents
        :type data : str
        :param options: ignored
        :type options : dict
        """

        geometries = []
        for line in data.split("\n"):
            extracted = self.line_to_geometry(line)
            if extracted:
                geometries.append(extracted)

        natoms = len(self.system.atoms)
        while geometries:
            g = geometries[:natoms]
            geometries = geometries[natoms:]
            self.geometry_history.append(g)

        self.geometry = self.geometry_history[-1]
        
    def run(self, host="localhost", options={}):
        """Run pdynamo using the the pdynamo-runner.py script on the given
        host. The pdynamo-runner.py gets geometry from an .xyz file and
        everything else from command line arguments. This runner needs to copy
        the geometry and the runner to the working directory before executing
        the runner.

        :param host: name of host where job should execute
        :type host : str
        :param options: ignored
        :type options : dict
        """

        run_params = self.get_run_config(host)
        workdir = self.backend + "-" + str(uuid.uuid1()).replace('-', '')[:16]
        path = "{0}/{1}/".format(self.tmpdir, workdir)
                
        #write .xyz geometry file to working directory with runner
        xyzdata = self.system.write("xyz")
        xyzfile = "{0}.xyz".format(workdir)
        xyzfull = "{0}{1}".format(path, xyzfile)
        self.write_file(xyzdata, xyzfull, host)

        #add geometry file specification and run location to command line
        log_file = path + xyzfile.replace(".xyz", ".log")
        args = {"path" : path, "input" : xyzfull, "output" : log_file}

        shell_cli = "python pdynamo-runner.py --xyzfile={input} --outfile={output} " + self.extras["cli_args"]
        shell_cli = shell_cli.format(**args)

        #copy runner to working directory
        here = os.path.dirname(os.path.dirname(__file__))
        src = "{0}/pdynamo-runner.py".format(here)
        with open(src) as runnerfile:
            runsource = runnerfile.read()
        dst = path + "pdynamo-runner.py"

        shell_script = path + "runpd.sh"
        self.write_file(shell_cli, shell_script, host)
        self.write_file(runsource, dst, host)

        cmd = self.deck.format(path=path)

        stdout, returncode = self.execute(cmd, host, bash_shell=True)
        
        self.stdout = stdout

        self.logdata = self.read_file(log_file, host)
        self.extract_last_energy(self.logdata)
        self.extract_heat_of_formation(self.logdata)
        self.extract_geometry(self.logdata)

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

        :param system: molecular system for energy calculation
        :type system : geoprep.System
        :param method: calculation method
        :type method : str
        :param options: additional keyword based control options
        :type options : dict
        :return: a pDynamo single point energy calculation job
        :rtype : cpinterface.Job
        """

        if method.startswith("semiempirical"):
            return self.make_semiempirical_job(system, method, "ENERGY",
                                               options=options)
        else:
            raise ValueError("pDynamo 7 does not support {0}".format(method))

    def make_opt_job(self, system, method, options={}):
        """Create an input specification for a geometry optimization
        calculation. Optimization goal may be to find minimum geometry or
        to find a saddle point.

        N.B.: There are actually many underlying geometry minimizers that are
        not yet exposed here. In addition to conjugate gradient minimization
        pDynamo has FIRE minimization, LBFGS minimization, quasi-Newton
        minimization, and steepest descent minimization.

        options:
         goal: minimize or saddle

        :param system: molecular system for energy calculation
        :type system : geoprep.System
        :param method: calculation method
        :type method : str
        :param options: additional keyword based control options
        :type options : dict
        :return: a pDynamo input for geometry optimization calculation
        :rtype : str
        """

        system = self.fragment_to_system(system)

        if method.startswith("semiempirical"):
            return self.make_semiempirical_job(system, method, "OPT",
                                               options=options)
        else:
            raise ValueError("pDynamo does not support {0}".format(method))

    def make_semiempirical_job(self, system, method, runtyp,
                               options={}):
        """Create a semiempirical input specification for a calculation.
        pDynamo supports AM1, MNDO, PM3, PM6, RM1, AM1/d-PhoT, PDDG/PM3,
        PDDG/MNDO

        :param system: molecular system for calculation
        :type system : geoprep.System
        :param method: a semiempirical calculation method
        :type method : str
        :param options: additional keyword based control options
        :type options : dict
        :return: a PDynamo semiempirical job
        :rtype : Job
        """

        defaults = {"reference" : "rhf", "goal" : "minimize"}
        options = dict(defaults.items() + options.items())
        goal = options.get("goal")

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

        elif runtyp == "OPT":
            controls.append("--jobtype=geometry")
            controls.append("--goal={0}".format(goal))

        controls.append("--charge={0}".format(system.charge))

        #The "deck" here is just a series of command line arguments for
        #pdynamo-runner.py. In the actual job runner we will add the name
        #of the geometry file, create a shell script, and copy the runner.
        args = " ".join([c for c in controls if c])
        deck = "cd {path} && bash ./runpd.sh "

        job = PDynamoJob(deck=deck, system=system, extras={"cli_args" : args})

        return job

