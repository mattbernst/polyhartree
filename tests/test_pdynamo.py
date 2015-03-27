#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_pdynamo
    ~~~~~~~~~~~~~~

    Test pDynamo implementations for semiempirical energy, minimum geometry,
    transition state geometry, and frequencies.
"""

import sys
import unittest
import geoprep
import pdynamo

class PDTestCase(unittest.TestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = pdynamo.PDynamo()

    def tearDown(self):
        pass

    def test_energy_pm3_methylium(self):
        methylium = self.G.make_system("[CH3+]")
        job = self.C.make_energy_job(methylium, "semiempirical:pm3")
        job.run_local()
        self.assertAlmostEqual(-5.641487, job.energy, places=5)
        self.assertAlmostEqual(0.408868, job.heat_of_formation, places=5)

    def test_energy_pm3_carbanide(self):
        carbanide = self.G.make_system("[CH3-]")
        job = self.C.make_energy_job(carbanide, "semiempirical:pm3")
        job.run_local()
        self.assertAlmostEqual(-5.967391, job.energy, places=5)
        self.assertAlmostEqual(0.082962, job.heat_of_formation, places=5)

    def test_energy_pm3_methyl_radical(self):
        mradical = self.G.make_system("[CH3]")
        job = self.C.make_energy_job(mradical, "semiempirical:pm3")
        self.assertEqual("Forcing UHF for multiplicity 2", self.C.messages[0])
        job.run_local()
        self.assertAlmostEqual(-6.005549, job.energy, places=5)
        self.assertAlmostEqual(0.044801, job.heat_of_formation, places=5)

    def test_energy_pm3_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pm3")
        job.run_local()
        self.assertAlmostEqual(-6.634475, job.energy, places=5)
        self.assertAlmostEqual(-0.020660, job.heat_of_formation, places=5)

    def test_energy_mndo_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:mndo")
        job.run_local()
        self.assertAlmostEqual(-6.801557, job.energy, places=5)
        self.assertAlmostEqual(-0.018602, job.heat_of_formation, places=5)

    def test_energy_am1_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:am1")
        job.run_local()
        self.assertAlmostEqual(-6.732507, job.energy, places=5)
        self.assertAlmostEqual(-0.012914, job.heat_of_formation, places=5)

    def test_energy_pm6_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pm6")
        job.run_local()
        self.assertAlmostEqual(-6.510840, job.energy, places=5)
        self.assertAlmostEqual(-0.019537, job.heat_of_formation, places=5)

    def test_energy_rm1_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:rm1")
        job.run_local()
        self.assertAlmostEqual(-6.716257, job.energy, places=5)
        self.assertAlmostEqual(-0.022076, job.heat_of_formation, places=5)

    def test_energy_pddgmndo_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pddg/mndo")
        job.run_local()
        self.assertAlmostEqual(-6.949286, job.energy, places=5)
        self.assertAlmostEqual(-0.026589, job.heat_of_formation, places=5)

    def test_energy_pddgpm3_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pddg/pm3")
        job.run_local()
        self.assertAlmostEqual(-6.727189, job.energy, places=5)
        self.assertAlmostEqual(-0.025640, job.heat_of_formation, places=5)

    def test_energy_am1dphot_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:am1/d-phot")
        job.run_local()
        self.assertAlmostEqual(-6.654096, job.energy, places=5)
        self.assertAlmostEqual(-0.002390, job.heat_of_formation, places=5)

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
        result = runSuite(PDTestCase, name = test_name)

    else:
        result = runSuite(PDTestCase)

    return result

if __name__ == "__main__":
    runTests()
