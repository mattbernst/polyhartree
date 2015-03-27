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
import gamess_us

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

    def test_bad_input_error(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pm3")

        #introduce an error in the input deck: misspell SCFTYP as SCFTYPE
        job.deck = job.deck.replace("SCFTYP", "SCFTYPE")

        self.assertEqual("begin", job.runstate)
        job.run_local()

        self.assertEqual("error", job.runstate)


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
