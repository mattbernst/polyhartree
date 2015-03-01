#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_mopac
    ~~~~~~~~~~~~~~

    Test MOPAC 7 implementations for energy, minimum geometry,
    transition state geometry, and frequencies.
"""

import sys
import unittest
from cinfony import pybel
import geoprep
import mopac7

class MOPACTestCase(unittest.TestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = mopac7.Mopac7()

    def tearDown(self):
        pass

    def test_methylium_energy_pm3(self):
        methylium = self.G.make_mol("[CH3+]")
        job = self.C.make_energy_job(methylium, "semiempirical:pm3")
        job.run_local()
        self.assertAlmostEqual(-10.85734, job.energy, places=4)
        self.assertAlmostEqual(0.408868, job.heat_of_formation, places=5)

    def test_carbanide_energy_pm3(self):
        carbanide = self.G.make_mol("[CH3-]")
        job = self.C.make_energy_job(carbanide, "semiempirical:pm3")
        job.run_local()
        self.assertAlmostEqual(-11.18324, job.energy, places=4)
        self.assertAlmostEqual(0.082962, job.heat_of_formation, places=5)

    def test_methyl_radical_energy_pm3(self):
        mradical = self.G.make_mol("[CH3]")
        job = self.C.make_energy_job(mradical, "semiempirical:pm3")
        job.run_local()
        self.assertAlmostEqual(-11.21848, job.energy, places=4)
        self.assertAlmostEqual(0.047725, job.heat_of_formation, places=5)

    def test_methane_energy_pm3(self):
        methane = self.G.make_mol("C")
        job = self.C.make_energy_job(methane, "semiempirical:pm3")
        job.run_local()
        self.assertAlmostEqual(-14.11915, job.energy, places=4)
        self.assertAlmostEqual(-0.020660, job.heat_of_formation, places=5)

    def test_methane_energy_mndo(self):
        methane = self.G.make_mol("C")
        job = self.C.make_energy_job(methane, "semiempirical:mndo")
        job.run_local()
        self.assertAlmostEqual(-14.34421, job.energy, places=4)
        self.assertAlmostEqual(-0.018578, job.heat_of_formation, places=5)

    def test_methane_energy_am1(self):
        methane = self.G.make_mol("C")
        job = self.C.make_energy_job(methane, "semiempirical:am1")
        job.run_local()
        self.assertAlmostEqual(-14.19059, job.energy, places=4)
        self.assertAlmostEqual(-0.012894, job.heat_of_formation, places=5)

    def test_methane_energy_mindo3(self):
        methane = self.G.make_mol("C")
        job = self.C.make_energy_job(methane, "semiempirical:mindo/3")
        job.run_local()
        self.assertAlmostEqual(-14.04530, job.energy, places=4)
        self.assertAlmostEqual(-0.009680, job.heat_of_formation, places=5)

    def test_error(self):
        #This doesn't work. At least the job status says so.
        #"THE FIRST THREE ATOMS MUST NOT LIE IN A STRAIGHT LINE"
        P2 = self.G.make_mol("P#P")
        job = self.C.make_energy_job(P2, "semiempirical:pm3")
        job.run_local()
        self.assertEqual("error", job.runstate)

    def test_element_error(self):
        #Try to create a system containing an element unparameterized for
        #MINDO/3. Should raise an error.
        #N.B.: The generated input file would actually run in practice,
        #rather contrary to the Mopac 7 manual's parameterization information,
        #but it seems best to err on the side of caution.
        SbCl3 = self.G.make_mol("Cl[Sb](Cl)Cl")
        self.assertRaises(ValueError, self.C.make_energy_job, SbCl3,
                          "semiempirical:mindo/3")

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
        test_name = sys.argv[2]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(MOPACTestCase, name = test_name)

    else:
        result = runSuite(MOPACTestCase)

    return result

if __name__ == "__main__":
    runTests()
