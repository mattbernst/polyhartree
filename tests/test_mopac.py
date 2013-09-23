#!/usr/bin/env python
# -*- mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
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
import gamess_us

class MOPACTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_methylium(self):
        G = geoprep.Geotool()
        C = gamess_us.GAMESSUS()
        methylium = G.make_mol("[CH3+]")
        job = C.make_energy_job(methylium, 'semiempirical:pm3')

    def test_carbanide(self):
        G = geoprep.Geotool()
        C = gamess_us.GAMESSUS()
        carbanide = G.make_mol("[CH3-]")
        job = C.make_energy_job(carbanide, 'semiempirical:pm3')

    def test_methyl_radical(self):
        G = geoprep.Geotool()
        C = gamess_us.GAMESSUS()
        mradical = G.make_mol("[CH3]")
        job = C.make_energy_job(mradical, 'semiempirical:pm3')

    def test_ethane(self):
        G = geoprep.Geotool()
        C = gamess_us.GAMESSUS()
        ethane = G.make_mol("CC")
        job = C.make_energy_job(ethane, 'semiempirical:pm3')
        import pdb; pdb.set_trace()


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
