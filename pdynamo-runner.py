# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import sys
import argparse
import os
import StringIO
import yaml
try:
    from pBabel import XYZFile_ToSystem
    from pCore import TextLogFileWriter
    from pMolecule import DIISSCFConverger, ElectronicState, QCModelMNDO
except ImportError:
    msg = "Can't find pDynamo components. Ensure that it is installed and on your PYTHON_PATH. http://www.pdynamo.org/\n"
    sys.stderr.write(msg)
    sys.exit(1)

controls = [("--xyzfile", {"help" : "XYZ file containing system atoms and geometry"}),
            ("--jobtype", {"help" : "Kind of job to run: energy"}),
            ("--method", {"help" : "Method for calculation: semiempirical:am1, semiempirical:mndo, semiempirical:pm3, semiempirical:pm6, semiempirical:rm1"}),
            ("--charge", {"help" : "System charge", "default" : 0, "type" : int}),
            ("--multiplicity", {"help" : "System multiplicity", "default" : 1, "type" : int}),
            ("--restricted", {"help" : "Restricted or unrestricted Hartree-Fock calculation (RHF or UHF)", "default" : "RHF"})]

class PDRunner(object):
    """A wrapper intended to call pDynamo for common computational chemistry
    tasks similar to the way one can call e.g. GAMESS or NWChem from the
    command line.

    Since pDynamo does not have its own job description language, directives
    are passed via command line arguments.
    """
    
    def __init__(self, *args, **kw):
        """Initialize runner with arguments from command line. Execute a job
        based on the arguments.

        @param *args: currently unused
        @type *args : tuple
        @param **kw: all keywords controlling job execution
        @type **kw : dict
        """

        self.logfile = StringIO.StringIO()
        self.logger = TextLogFileWriter()
        self.logger.file = self.logfile
        
        system = XYZFile_ToSystem(kw["xyzfile"])

        jobtype = kw["jobtype"]
        charge = kw["charge"]
        multiplicity = kw["multiplicity"]
        #uhf means unrestricted calculation, rhf means restricted
        restricted = {"uhf" : False, "rhf" : True}[kw["restricted"].lower()]

        method = kw["method"]
        if method.startswith("semiempirical:"):
            method_name = method.split(":")[1]
        else:
            raise ValueError("Unrecognized method {0}".format(method))

        if jobtype == "energy":
            system = self.run_spe_semiempirical(system, method_name, charge,
                                                multiplicity, restricted)
            hof = self.get_heat_of_formation(system, method_name)
            t = self.logger.GetTable()
            t.Start()
            v = "Heat of Formation: {:.6f} kcal/mol".format(hof)
            t.Title(v)
            t.Stop()
            sys.stdout.write(self.logfile.getvalue())
        else:
            raise ValueError("Unrecognized job type {0}".format(jobtype))

    def run_spe_semiempirical(self, system, method_name, charge, multiplicity,
                              restricted):
        """Run a single point energy calculation for a semiempirical method
        like AM1 or PM6.

        @param system: input system with geometry
        @type system : pMolecule.System.System
        @param method_name: name of semiempirical model to apply
        @type method_name : str
        @param charge: system charge
        @type charge : int
        @param multiplicity: system multiplicity
        @type multiplicity : int
        @param restricted: whether to run a spin-restricted calculation
        @type restricted : bool
        """

        es = ElectronicState(charge=charge, multiplicity=multiplicity)
        converger = DIISSCFConverger(densityTolerance=1.0e-6,
                                     maximumSCFCycles=999)

        qcmodel = QCModelMNDO(method_name, converger=converger,
                              isSpinRestricted=restricted)

        system.electronicState = es
        system.DefineQCModel(qcmodel)
        system.Energy(log=self.logger)
        
        return system

    def get_se_atom_isol(self, semiempirical_method, element):
        """Get E(isol-atom) for the given semiempirical method and element.
        The value is returned in eV.

        @param semiempirical_method: a method name, like "am1"
        @type semiempirical_method : str
        @param element: an element symbol, like "Cl"
        @type element : str
        @return: energy
        @rtype : float
        """
        
        efile = "{0}.yaml".format(element)
        parts = [os.environ["PDYNAMO_PARAMETERS"], "mndoParameters",
                 semiempirical_method, efile]

        f = os.path.sep.join(parts)
        with open(f) as infile:
            data = yaml.load(infile)

        params = data["Scalar Parameters"]

        #need to get per-atom E-isol and E-atom,
        #but E-isol is eV and E-atom is kcal/mol
        eisol = [r[1] for r in params if r[0] == "eisol"][0]
        eatom = [r[1] for r in params if r[0] == "eheat"][0]
        eatom_ev = eatom * 0.04336410268575565

        return (eisol - eatom_ev)

    def get_heat_of_formation(self, converged_system, semiempirical_method):
        """Calculate the heat of formation for a system that has reached
        electronic convergence, based on a semiempirical model.

        The heat of formation is calculated as explained at
        http://openmopac.net/manual/SCF_calc_hof.html

        @param converged_system: a system with energy already calculated
        @type converged_system : pMolecule.System.System
        @param semiempirical_method: name of SE method for parameter lookup
        @type semiempirical_method : str
        @return: heat of formation in kcal/mol
        @rtype : float
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

        history = self.logfile.getvalue()
        for line in history.split("\n"):
            if "Electronic Energy" in line and "Nuclear Energy" in line:
                pieces = line.split()
                electronic = float(pieces[3])
                nuclear = float(pieces[7])

        atnos = [a.atomicNumber for a in converged_system.atoms.items]
        elems = [elements[a - 1] for a in atnos]
        aes = [self.get_se_atom_isol(semiempirical_method, e) for e in elems]

        #For reasons unknown the electronic and nuclear energies are somewhat
        #different vs Mopac7, e.g. in a methane calculation with identical
        #geometries Mopac7 gives -387.187423 eV electronic 206.676675 eV nuclear
        #and this gives -389.187352 eV electronic 206.098290 eV nuclear
        electronic_ev = electronic * 27.21138505
        nuclear_ev = nuclear * 27.21138505
        total = (electronic + nuclear) - (sum(aes) / 27.21138505)
        total_kcalm = total * 627.509469

        return total_kcalm

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    keys = []
    for c in controls:
        #get name without -- prefix, e.g. "method" instead of "--method"
        keys.append(c[0][2:])
        parser.add_argument(c[0], **c[1])
        
    args = parser.parse_args()
    settings = {}
    for key in keys:
        settings[key] = getattr(args, key)

    pd = PDRunner(**settings)
