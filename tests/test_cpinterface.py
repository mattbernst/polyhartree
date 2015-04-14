#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_cpinterface
    ~~~~~~~~~~~~~~

    Test chemical program interface code that is not tied to any one specific
    back-end.
"""

import sys
import unittest
import cpinterface
import geoprep

class CPITestCase(unittest.TestCase):
    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = cpinterface.MolecularCalculator()

    def tearDown(self):
        pass

    def test_missing_basis_name(self):
        #error given for atoms lacking a basis set assignment
        methylium = self.G.make_fragment("[CH3+]")
        methane = self.G.make_fragment("C")
        methylium.set_basis_name("cc-pVDZ")

        s = geoprep.System([methylium, methane])
        self.assertRaises(ValueError, self.C.get_basis_data, s,
                          {"basis_format" : "gamess-us"})

    def test_bad_basis_name(self):
        #error given for atoms with unknown basis set name
        methylium = self.G.make_fragment("[CH3+]")
        methane = self.G.make_fragment("C")
        methylium.set_basis_name("cc-pVDZ")
        methane.set_basis_name("bozo-basis")

        s = geoprep.System([methylium, methane])
        self.assertRaises(ValueError, self.C.get_basis_data, s,
                          {"basis_format" : "gamess-us"})

    def test_mixed_basis_representation(self):
        #error given for attempting to mix spherical and cartesian basis
        #functions in a single system
        methylium = self.G.make_fragment("[CH3+]")
        methane = self.G.make_fragment("C")
        methylium.set_basis_name("cc-pVDZ")
        methane.set_basis_name("6-31G")

        s = geoprep.System([methylium, methane])
        self.assertRaises(ValueError, self.C.get_basis_data, s,
                          {"basis_format" : "gamess-us"})

    def test_missing_basis_elements(self):
        #error given for elements unparameterized by chosen basis set
        lewisite = self.G.make_fragment("Cl[As](Cl)\C=C\Cl")
        lewisite.set_basis_name("TZ (Dunning)")
        
        s = geoprep.System([lewisite])
        self.assertRaises(ValueError, self.C.get_basis_data, s,
                          {"basis_format" : "gamess-us"})

    def test_basis_retrieval(self):
        #basic error-free basis retrieval
        lewisite = self.G.make_fragment("Cl[As](Cl)\C=C\Cl")
        lewisite.set_basis_name("cc-pVDZ")
        chlorines = lewisite.select("[Cl]", hydrogen="exclude")
        
        lewisite.set_basis_name("cc-pVTZ", chlorines)
        
        s = geoprep.System([lewisite])
        bd = self.C.get_basis_data(s, {"basis_format" : "gamess-us"})
        self.assertEqual("spherical", bd["spherical_or_cartesian"])
        
        data = bd["data"]
        self.assertEqual(["cc-pVDZ", "cc-pVTZ"], sorted(data.keys()))

        tzd = "".join(data["cc-pVTZ"].values())
        self.assertTrue("CHLORINE" in tzd)

        dzd = "".join(data["cc-pVDZ"].values())
        for name in ["HYDROGEN", "CARBON", "ARSENIC"]:
            self.assertTrue(name in dzd)

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
        result = runSuite(CPITestCase, name = test_name)

    else:
        result = runSuite(CPITestCase)

    return result

if __name__ == '__main__':
    runTests()
