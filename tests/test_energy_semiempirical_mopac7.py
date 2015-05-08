#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_energy_semiempirical_mopac7
    ~~~~~~~~~~~~~~

    Test MOPAC 7 semiempirical implementation for energy.
"""

import sys
import unittest
import geoprep
from adapters import mopac7
from tests import energy_semiempirical as es
from tests import reference_values

class MOPACEnergyTestCase(es.SemiempiricalEnergyTestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = mopac7.Mopac7()

    def tearDown(self):
        pass

    def test_energy_mindo3_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:mindo/3")
        job.run()
        self.assertAlmostEqual(reference_values.methane_mindo3_hof,
                               job.heat_of_formation, places=5)

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
        result = runSuite(MOPACEnergyTestCase, name = test_name)

    else:
        result = runSuite(MOPACEnergyTestCase)

    return result

if __name__ == "__main__":
    runTests()
