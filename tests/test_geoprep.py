# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_geoprep
    ~~~~~~~~~~~~~~

    Test backend-independent functionality: fragment, system, basis name
    assignment...
"""

import sys
import unittest
from cinfony import pybel
import geoprep

class GTestCase(unittest.TestCase):

    def setUp(self):
        self.G = geoprep.Geotool()

    def tearDown(self):
        pass

    def test_spin_assignment(self):
        #create triplet oxygen
        oxygen = self.G.make_fragment("[O][O]")
        self.assertEqual(3, oxygen.spin)
        
        #make it singlet
        oxygen.spin = 1
        self.assertEqual(1, oxygen.spin)

    def test_nelec(self):
        #count the number of electrons

        methyl_radical = self.G.make_fragment("[CH3]")
        methylium = self.G.make_fragment("[CH3+]")
        carbanide = self.G.make_fragment("[CH3-]")

        self.assertEqual(9, methyl_radical.nelec)
        self.assertEqual(methyl_radical.nelec + 1, carbanide.nelec)
        self.assertEqual(methyl_radical.nelec - 1, methylium.nelec)


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
        result = runSuite(GTestCase, name = test_name)

    else:
        result = runSuite(GTestCase)

    return result

if __name__ == "__main__":
    runTests()
