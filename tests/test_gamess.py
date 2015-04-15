#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_gamess
    ~~~~~~~~~~~~~~

    Test GAMESS implementations for energy, minimum geometry,
    transition state geometry, and frequencies.
"""

import sys
import unittest
import geoprep
from adapters import gamess_us

class GAMESSTestCase(unittest.TestCase):
    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = gamess_us.GAMESSUS()

    def tearDown(self):
        pass

    def test_energy_rohf_uhf_pm3(self):
        #compare UHF and ROHF across different radicals for heat of formation
        
        smiles = {"methyl_radical" : "[CH3]",
                  "phenyl_radical" : "c1cc[c]cc1"}
        jobs = {}

        for key in smiles:
            s = smiles[key]
            mol = self.G.make_system(s)
            for reference in ["uhf", "rohf"]:
                j = self.C.make_energy_job(mol, "semiempirical:pm3",
                                           options={"reference" : reference})
                j.run_local()
                name = "{0}-{1}".format(key, reference)
                jobs[name] = j

        hofs = dict([(j, jobs[j].heat_of_formation) for j in jobs])

        self.assertAlmostEqual(0.126826, hofs["phenyl_radical-rohf"], places=5)
        self.assertAlmostEqual(0.119951, hofs["phenyl_radical-uhf"], places=5)
        self.assertAlmostEqual(0.047725, hofs["methyl_radical-rohf"], places=5)
        self.assertAlmostEqual(0.044800, hofs["methyl_radical-uhf"], places=5)

    def test_energy_pm3_methylium(self):
        methylium = self.G.make_system("[CH3+]")
        job = self.C.make_energy_job(methylium, "semiempirical:pm3")
        job.run_local()
        self.assertAlmostEqual(-5.641425, job.energy, places=5)
        self.assertAlmostEqual(0.408868, job.heat_of_formation, places=5)

    def test_energy_pm3_carbanide(self):
        carbanide = self.G.make_system("[CH3-]")
        job = self.C.make_energy_job(carbanide, "semiempirical:pm3")
        job.run_local()
        self.assertAlmostEqual(-5.967321, job.energy, places=5)
        self.assertAlmostEqual(0.082962, job.heat_of_formation, places=5)

    def test_energy_pm3_methyl_radical(self):
        methyl_radical = self.G.make_system("[CH3]")
        job = self.C.make_energy_job(methyl_radical, "semiempirical:pm3")
        self.assertEqual("Forcing UHF for multiplicity 2", self.C.messages[0])
        job.run_local()
        self.assertAlmostEqual(-6.005483, job.energy, places=5)
        self.assertAlmostEqual(0.044800, job.heat_of_formation, places=5)

    def test_energy_pm3_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, 'semiempirical:pm3')
        job.run_local()
        self.assertAlmostEqual(-6.634400, job.energy, places=5)
        self.assertAlmostEqual(-0.020660, job.heat_of_formation, places=5)

    def test_energy_mndo_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:mndo")
        job.run_local()
        self.assertAlmostEqual(-6.801455, job.energy, places=5)
        self.assertAlmostEqual(-0.018578, job.heat_of_formation, places=5)

    def test_energy_am1_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:am1")
        job.run_local()
        self.assertAlmostEqual(-6.732408, job.energy, places=5)
        self.assertAlmostEqual(-0.012894, job.heat_of_formation, places=5)

    def test_energy_rm1_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:rm1")
        job.run_local()
        self.assertAlmostEqual(-6.716163, job.energy, places=5)
        self.assertAlmostEqual(-0.022059, job.heat_of_formation, places=5)

    def test_energy_scf_methane(self):
        #very basic minimal basis set test for methane
        methane = self.G.make_fragment("C")
        methane.set_basis_name("3-21G")
        job = self.C.make_energy_job(methane, "hf:rhf")
        job.run_local()
        self.assertAlmostEqual(-39.976642, job.energy, places=5)

    def test_energy_rohf_uhf_scf_methane(self):
        #compare UHF and ROHF across methyl radicals for HF energy
        
        smiles = {"methyl_radical" : "[CH3]"}
        jobs = {}

        for key in smiles:
            s = smiles[key]
            fragment = self.G.make_fragment(s)
            fragment.set_basis_name("cc-pVDZ")

            for reference in ["uhf", "rohf"]:
                j = self.C.make_energy_job(fragment,
                                           "hf:{0}".format(reference))
                j.run_local()
                name = "{0}-{1}".format(key, reference)
                jobs[name] = j

        energies = dict([(j, jobs[j].energy) for j in jobs])

        self.assertAlmostEqual(-39.5638133, energies["methyl_radical-rohf"],
                               places=6)
        self.assertAlmostEqual(-39.5638132, energies["methyl_radical-uhf"],
                               places=6)

    def test_bad_input_error(self):
        methane = self.G.make_fragment("C")
        job = self.C.make_energy_job(methane, "semiempirical:pm3")

        #introduce an error in the input deck: misspell SCFTYP as SCFTYPE
        job.deck = job.deck.replace("SCFTYP", "SCFTYPE")

        self.assertEqual("begin", job.runstate)
        job.run_local()

        self.assertEqual("error", job.runstate)

    def test_reformat_long_line(self):
        #verify splitting of long directives so as to comfortably fit within
        #the 80 character per line limit of GAMESS-US
        maxlen = 50
        line = " $BASIS basnam(1)=ClBAS0, AsBAS0, ClBAS0, CBAS0, CBAS0, ClBAS0, HBAS0, HBAS0 $END"
        pieces, begin, end = self.C.reformat_long_line(line, " $BASIS", "$END",
                                                       maxlen=maxlen)
        for piece in pieces:
            self.assertTrue(len(piece) <= maxlen)

    def test_prepare_basis_data(self):
        #test generation of inline basis set data with water
        #should generate one basis set assignment for oxygen and two for
        #hydrogen, since one hydrogen gets assigned a different basis
        expected_comments = "!Basis set assignments:\n!\tO 3-21G\n!\tH 6-31G\n!\tH 3-21G"
        expected_labels = ["$HBAS0", "$HBAS1", "$OBAS0"]

        water = self.G.make_fragment("O")
        water.set_basis_name("3-21G")
        hydrogens = water.select("[O]", hydrogen="only")
        water.set_basis_name("6-31G", hydrogens[:1])
        
        s = geoprep.System([water])
        job = self.C.make_energy_job(s, "hf:rhf")
        self.assertTrue(expected_comments in job.deck)
        for label in expected_labels:
            self.assertTrue(label in job.deck)


def runSuite(cls, verbosity=2, name=None):
    """Run a unit test suite and return status code.

    @param cls: class that the suite should be constructed from
    @type cls : class
    @param verbosity: verbosity level to pass to test runner
    @type verbosity : int
    @param name: name of a specific test in the suite to run
    @type name : str
    @return: unit test run status code
    @rtype : int
    """
    try: 
        if name:
            suite = unittest.makeSuite(cls, name)
        else:
            suite = unittest.makeSuite(cls)
            
        return unittest.TextTestRunner(verbosity=verbosity).run(suite)
    
    except SystemExit:
        pass

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(GAMESSTestCase, name = test_name)

    else:
        result = runSuite(GAMESSTestCase)

    return result

if __name__ == '__main__':
    runTests()
