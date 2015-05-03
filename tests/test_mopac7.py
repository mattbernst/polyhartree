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
import geoprep
from adapters import mopac7
from tests import reference_values

class MOPACTestCase(unittest.TestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = mopac7.Mopac7()

    def tearDown(self):
        pass

    def test_energy_pm3_methylium(self):
        methylium = self.G.make_system("[CH3+]")
        job = self.C.make_energy_job(methylium, "semiempirical:pm3")
        job.run()
        self.assertAlmostEqual(reference_values.methylium_pm3_hof,
                               job.heat_of_formation, places=5)

    def test_energy_pm3_carbanide(self):
        carbanide = self.G.make_system("[CH3-]")
        job = self.C.make_energy_job(carbanide, "semiempirical:pm3")
        job.run()
        self.assertAlmostEqual(reference_values.carbanide_pm3_hof,
                               job.heat_of_formation, places=5)

    def test_energy_pm3_methyl_radical(self):
        mradical = self.G.make_system("[CH3]")
        job = self.C.make_energy_job(mradical, "semiempirical:pm3")
        self.assertEqual("Forcing UHF for multiplicity 2", self.C.messages[0])
        job.run()
        self.assertAlmostEqual(reference_values.methyl_uhf_pm3_hof,
                               job.heat_of_formation, places=5)

    def test_energy_pm3_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pm3")
        job.run()
        self.assertAlmostEqual(reference_values.methane_pm3_hof,
                               job.heat_of_formation, places=5)

    def test_energy_mndo_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:mndo")
        job.run()
        self.assertAlmostEqual(reference_values.methane_mndo_hof,
                               job.heat_of_formation, places=5)

    def test_energy_am1_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:am1")
        job.run()
        self.assertAlmostEqual(reference_values.methane_am1_hof,
                               job.heat_of_formation, places=5)

    def test_energy_mindo3_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:mindo/3")
        job.run()
        self.assertAlmostEqual(reference_values.methane_mindo3_hof,
                               job.heat_of_formation, places=5)

    def test_bad_input_error(self):
        #This doesn't work. At least the job status says so.
        #"THE FIRST THREE ATOMS MUST NOT LIE IN A STRAIGHT LINE"
        P2 = self.G.make_system("P#P")
        job = self.C.make_energy_job(P2, "semiempirical:pm3")
        job.run()
        self.assertEqual("error", job.runstate)

    def test_element_error(self):
        #Try to create a system containing an element unparameterized for
        #MINDO/3. Should raise an error.
        #N.B.: The generated input file would actually run in practice,
        #rather contrary to the Mopac 7 manual's parameterization information,
        #but it seems best to err on the side of caution.
        stibine = self.G.make_system("[SbH3]")
        self.assertRaises(ValueError, self.C.make_energy_job, stibine,
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
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(MOPACTestCase, name = test_name)

    else:
        result = runSuite(MOPACTestCase)

    return result

if __name__ == "__main__":
    runTests()
