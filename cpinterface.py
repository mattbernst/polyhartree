# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import os
import cPickle as pickle
import sys
import yaml
import sharedutilities
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

class Job(sharedutilities.Utility, Messages):
    def __init__(self, deck="", system=None, runstate="begin", tmpdir="/tmp",
                 extras={}):
        #states: begin, running, complete, error
        self.runstate = runstate
        self.system = system
        self.deck = deck
        self.stdout = ""
        self.logdata = ""
        self.tmpdir = tmpdir
        self.energy = None
        self.heat_of_formation = None
        self.geometry = []
        self.geometry_history = []
        self.messages = []
        self.extras = extras
        #use location of script to find location of configs 
        self.here = os.path.dirname(__file__)
        #2 days
        self.timeout = 86400 * 2

    def get_run_config(self, host):
        """Get the run configuration for the backend used by this job.
        If config/runners.yaml is present it will take precedence over
        the default config/runners-default.yaml.

        :param host: name of host to use for job execution
        :type host : str
        :return: backend run configuration
        :rtype : dict
        """

        data = None
        
        for fname in ["config/runners.yaml", "config/runners-default.yaml"]:
            try:
                with open(fname) as infile:
                    data = yaml.safe_load(infile)
                    break
            except IOError:
                continue

        if data is None:
            raise IOError("Unable to locate config file runners.yaml or runners-default.yaml")
        configs = data[host]
        config = None
        for c in configs:
            eflag = c["enabled"].lower()
            if c["program"] == self.backend and eflag in ("y", "true"):
                config = c
                break

        if config is None:
            raise KeyError("Could not find enabled backend for {0} on {1}".format(self.backend, host))
            
        return config

    def run(self, host="localhost", options={}):
        raise NotImplementedError

    def ansible_run(self, module_name, module_args, host, complex_args={}):
        """Synchronously run an ansible command on a single host.

        N.B.: You may need to edit ~/.bashrc on the machines where you are
        running chemistry programs, if your .bashrc contains important setup
        logic like extending PATH or setting other environment varables.
        On most distributions the .bashrc is largely disabled when a
        non-interactive shell is being executed. Ansible foolishly does not
        permit you to run shell commands in full interactive mode.

        Look in .bashrc for a line like this:
        # If not running interactively, don't do anything

        and then comment out the logic below it.

        :param module_name: the ansible module to run
        :type module_name : str
        :param module_args: the ansible module arguments to pass
        :type module_args : dict
        :param host: the host to target with ansible
        :type host : str
        :param complex_args: more optional arguments to pass to the runner
        :type complex_args : dict
        """

        self.load_ansible()

        loaded_host = self.inventory.get_host(host)
        if loaded_host is None:
            raise KeyError("Host {0} unknown -- check your ansible-hosts configuration".format(host))
        run_args = {"module_name" : module_name, "module_args" : module_args,
                    "forks" : 1, "subset" : host, "timeout" : self.timeout,
                    "inventory" : self.inventory,
                    "complex_args" : complex_args}

        r = ansible.runner.Runner(**run_args)
        result = r.run()
        self.log_ansible_errors(result)
        
        return result

    def log_ansible_errors(self, result):
        """Display and store any ansible command errors that might be
        returned.

        :param result: ansible output from calling runner.run()
        :type result : dict
        """

        for host in result["dark"]:
            msg = result["dark"][host].get("msg", "")
            message = "{0} : {1}".format(host, msg)
            sys.stderr.write(message + "\n")
            self.log(message)

    def write_file(self, data, filename, host):
        """Write the given data to filename on host. Also create the
        destination directory first if it does not yet exist. If host is
        anything other than localhost, use ansible for remote write.

        :param data: data to write
        :type data : str
        :param filename: name of file to write, with absolute path prepended
        :type filename : str
        """

        if host == "localhost":
            dirname = os.path.dirname(filename)
            if not os.path.exists(dirname):
                os.makedirs(dirname)

            with open(filename, "wb") as outfile:
                outfile.write(data)

        else:
            #copy departing file to temporary location instead of directly
            #passing data to ansible as 'content' parameter, because 'content'
            #tries to examine data and can get confused
            dirname = os.path.dirname(filename) + "/departing/"
            if not os.path.exists(dirname):
                os.makedirs(dirname)

            localname = dirname + os.path.basename(filename)
            with open(localname, "wb") as tempfile:
                tempfile.write(data)
                
            mkdir_args = {"path" : os.path.dirname(filename),
                          "state" : "directory"}
            r1 = self.ansible_run("file", mkdir_args, host)
            copy_args = {"src" : localname, "dest" : filename}
            r2 = self.ansible_run("copy", copy_args, host)

    def read_file(self, filename, host):
        """Read and return the data from filename on host. If the data is on
        a host other than localhost, copy it to the local machine first.

        :param filename: name of file to read, with absolute path prepended
        :type filename : str
        :return: file contents
        :rtype : str
        """

        if host == "localhost":
            with open(filename, "rb") as infile:
                data = infile.read()

        else:
            dirname = os.path.dirname(filename) + "/arriving/"
            if not os.path.exists(dirname):
                os.makedirs(dirname)

            destination = dirname + os.path.basename(filename)
            args = {"src" : filename, "dest" : destination, "flat" : True}
            self.ansible_run("fetch", args, host)
            try:
                with open(destination, "rb") as infile:
                    data = infile.read()
            except IOError:
                data = None

        return data

    def execute(self, cmd, host, stdin_data="", bash_shell=False, cwd=None):
        """Execute a command. Run locally if host is localhost or over ansible
        otherwise.

        :return: output, return code
        :rtype : tuple
        """

        if host == "localhost":
            output, rcode = self.execute_local(cmd, stdin_data=stdin_data,
                                               bash_shell=bash_shell, cwd=cwd)

        else:
            shell_binary = {"executable" : "/bin/bash"}
            result = self.ansible_run("shell", cmd, host,
                                      complex_args=shell_binary)
            run = result["contacted"].get(host, {})
            output = run.get("stdout")
            rcode = run.get("rc")

        return (output, rcode)

    def load_ansible(self):
        """Import ansible modules just before executing remote tasks. This
        way someone who just wants to run locally has one less dependency
        to install."""
        try:
            global ansible
            import ansible.inventory
            import ansible.runner
        except ImportError:
            msg = "You must have ansible installed to run jobs on remote machines. See http://www.ansible.com/home\n"
            sys.stderr.write(msg)
            sys.exit(1)

        try:
            self.inventory
        except AttributeError:
            
            for hostfile in ["config/ansible-hosts",
                             "config/ansible-hosts-default"]:
                full_file = "{0}/{1}".format(self.here, hostfile)
                if os.path.exists(full_file):
                    inv = ansible.inventory.Inventory(host_list=hostfile)
                    break

            self.inventory = inv

    def extract_last_energy(self, data, options={}):
        raise NotImplementedError

    def extract_heat_of_formation(self, data, options={}):
        raise NotImplementedError

    def match_line(self, line, pattern, indices):
        """Numericize a line and when it matches a type pattern, extract
        the whitespace separated components indexed by indices. This is
        useful for e.g. getting geometry out of log files.

        :param line: a line of data to process
        :type line : str
        :param pattern: types to match, e.g. [int, str, float, float, float]
        :type pattern : list
        :param indices: indices of components to extract and return from matches
        :type indices : list
        :return: matched components, if there was a match
        :rtype : list
        """

        matched = []
        n = self.numericize(line, force_float=False)
        if [type(k) for k in n] == pattern:
            matched = [n[k] for k in indices]

        return matched

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
        pattern, match_indices = self.geometry_matcher
        for line in data.split("\n"):
            extracted = self.match_line(line, pattern, match_indices)
            if extracted:
                geometries.append(extracted)

        elements = [sharedutilities.ELEMENTS[k.atomicnum - 1] for k in self.system.atoms]
        natoms = len(elements)

        while geometries:
            g = geometries[:natoms]
            geometries = geometries[natoms:]

            for j in range(natoms):
                g[j] = [elements[j]] + g[j]

            self.geometry_history.append(g)

        #remove any repeated geometries
        dupes = set()
        deduped = []
        for g in self.geometry_history:
            dumped = pickle.dumps(g)
            if dumped not in dupes:
                deduped.append(g)
                dupes.add(dumped)
        self.geometry_history = deduped

        self.geometry = self.geometry_history[-1]

class MolecularCalculator(Messages):
    def __init__(self, *args, **kw):
        self.messages = []
        
    def create_geometry(self, molecule, options={}):
        raise NotImplementedError

    def make_energy_job(self, system, method, options={}):
        raise NotImplementedError

    def fragment_to_system(self, sf):
        """It's convenient to be able to pass in a single-fragment system as
        a fragment without manually converting it to a system.

        :param sf: a system or a single fragment
        :type sf : geoprep.Fragment | geoprep.System
        :return: a system
        :rtype : geoprep.System
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

        :param system: molecular system
        :type system : geoprep.System
        :param options: ignored
        :type options : dict
        :return: basis data block
        :rtype : str
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

        :param reference: reference name to check
        :type reference : str
        :return: original reference name
        :rtype : str
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

        :param method: method name to check
        :type reference : str
        :return: original method name
        :rtype : str
        """
        
        if method in self.methods:
            return method
            
        else:
            raise ValueError("Unrecognized method {0}".format(repr(method)))

    def check_coordinates(self, coordinate_choice):
        """Check that coordinate system is supported by the quantum chemistry
        back-end in use.

        Some examples: cartesian, zmatrix

        :param coordinate_choice: coordinate system name to check
        :type reference : str
        :return: original coordinate system name
        :rtype : str
        """
        
        if coordinate_choice in self.coordinate_choices:
            return coordinate_choice

        else:
            raise ValueError("Unrecognized coordinates option {0}".format(repr(coordinate_choice)))
