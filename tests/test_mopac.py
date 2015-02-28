#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_mopac
    ~~~~~~~~~~~~~~

    Test MOPAC implementations for geometry optimization and energy.

"""

import sys
import unittest
from cinfony import pybel
import geoprep
import mopac7

class MOPACTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_methylium_energy(self):
        G = geoprep.Geotool()
        C = mopac7.Mopac7()
        methylium = G.make_mol("[CH3+]")
        job = C.make_energy_job(methylium, 'semiempirical:pm3')
        job.run_local()
        self.assertAlmostEqual(-10.85734, job.energy, places=4)
        self.assertAlmostEqual(0.408868, job.heat_of_formation,
                               places=5)

    def test_carbanide_energy(self):
        G = geoprep.Geotool()
        C = mopac7.Mopac7()
        carbanide = G.make_mol("[CH3-]")
        job = C.make_energy_job(carbanide, 'semiempirical:pm3')
        job.run_local()
        self.assertAlmostEqual(-11.18324, job.energy, places=4)
        self.assertAlmostEqual(0.082962, job.heat_of_formation, places=5)

    def test_methyl_radical_energy(self):
        G = geoprep.Geotool()
        C = mopac7.Mopac7()
        mradical = G.make_mol("[CH3]")
        job = C.make_energy_job(mradical, 'semiempirical:pm3')
        job.run_local()
        self.assertAlmostEqual(-11.22140, job.energy, places=4)
        self.assertAlmostEqual(0.044800, job.heat_of_formation, places=5)

    def test_methane_energy(self):
        G = geoprep.Geotool()
        C = mopac7.Mopac7()
        methane = G.make_mol("C")
        job = C.make_energy_job(methane, 'semiempirical:pm3')
        job.run_local()
        self.assertAlmostEqual(-14.11915, job.energy, places=4)
        self.assertAlmostEqual(-0.020660, job.heat_of_formation, places=5)

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

if __name__ == '__main__':
    runTests()
