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

    def test_spin_guess(self):
        #test spin multiplicity as guessed by underlying OpenBabel code
        methane = self.G.make_fragment("[CH4]")
        methyl_radical = self.G.make_fragment("[CH3]")
        methylene_diradical = self.G.make_fragment("[CH2]")

        self.assertEqual(1, methane.spin)
        self.assertEqual(2, methyl_radical.spin)

        #This could be singlet or triplet. OpenBabel assumes triplet. No way
        #to specify ambiguous multiplicity in standard SMILES.
        self.assertEqual(3, methylene_diradical.spin)

    def test_system_nelec(self):
        #count all electrons (sum of per-fragment nelec)
        methane = self.G.make_fragment("C")
        water = self.G.make_fragment("O")
        s = geoprep.System([methane, water])

        self.assertEqual(methane.nelec + water.nelec, s.nelec)

    def test_system_spin_guess(self):
        #spin multiplicity is inherited from the first fragment in the system
        #unless explicitly set
        methylene_diradical = self.G.make_fragment("[CH2]")
        water = self.G.make_fragment("O")
        s = geoprep.System([methylene_diradical, water])

        self.assertEqual(methylene_diradical.spin, s.spin)
        s.spin = 1
        self.assertEqual(1, s.spin)
        
    def test_system_charge(self):
        #system charge is the sum of fragment charges

        azanide = self.G.make_fragment("[NH2-]")
        methylium = self.G.make_fragment("[CH3+]")
        s = geoprep.System([azanide, methylium, azanide])

        self.assertEqual(-1, s.charge)

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
