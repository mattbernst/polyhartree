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
from adapters import pdynamo
from tests import reference_values

class PDTestCase(unittest.TestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = pdynamo.PDynamo()

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
                               job.heat_of_formation, places=4)

    def test_energy_am1_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:am1")
        job.run()
        self.assertAlmostEqual(reference_values.methane_am1_hof,
                               job.heat_of_formation, places=4)

    def test_energy_rm1_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:rm1")
        job.run()
        self.assertAlmostEqual(reference_values.methane_rm1_hof,
                               job.heat_of_formation, places=4)

    def test_energy_pm6_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pm6")
        job.run()
        self.assertAlmostEqual(reference_values.methane_pm6_hof,
                               job.heat_of_formation, places=5)

    def test_energy_pddgmndo_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pddg/mndo")
        job.run()
        self.assertAlmostEqual(reference_values.methane_pddgmndo_hof,
                               job.heat_of_formation, places=5)

    def test_energy_pddgpm3_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pddg/pm3")
        job.run()
        self.assertAlmostEqual(reference_values.methane_pddgpm3_hof,
                               job.heat_of_formation, places=5)

    def test_energy_am1dphot_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:am1/d-phot")
        job.run()
        self.assertAlmostEqual(reference_values.methane_am1dphot_hof,
                               job.heat_of_formation, places=5)

    def test_extract_geometry_unchanged(self):
        #extract geometry from an energy job, which should be basically the
        #same as the input
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pm3")
        job.run()

        reread = self.G.geolist_to_fragment(job.geometry)
        rmsd = self.G.align(methane.fragments[0], reread)["rmsd"]
        self.assertTrue(rmsd < 10**-5)


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
