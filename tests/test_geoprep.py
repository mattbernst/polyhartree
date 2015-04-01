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

    def test_system_elements(self):
        #unique elements in system, ordered by ascending atomic number
        
        halothane = self.G.make_fragment("C(C(F)(F)F)(Cl)Br")
        cysteine = self.G.make_fragment("C([C@@H](C(=O)O)N)S")
        s = geoprep.System([halothane, cysteine])
        expected = ["H", "C", "N", "O", "F", "S", "Cl", "Br"]
        
        self.assertEqual(expected, s.elements)

    def test_system_title(self):
        #system title can be set explicitly or inherits from first fragment
        name = "azanide"
        azanide = self.G.make_fragment("[NH2-]")
        azanide.title = name
        methylium = self.G.make_fragment("[CH3+]")
        s = geoprep.System([azanide, methylium, azanide])

        self.assertEqual(name, s.title)

        name2 = "amide ion"
        s.title = name2
        self.assertEqual(name2, s.title)

    def test_hydrogen_selectors(self):
        triethylamine = self.G.make_fragment("CCN(CC)CC")
        ethyl = "[#6][#6]"

        #heavy atoms (carbon) of ethyl groups alone
        e1 = [[0, 1], [3, 4], [5, 6]]
        #hydrogen atoms attached to ethyl groups, alone
        e2 = [[7, 8, 9, 10, 11], [12, 13, 14, 15, 16], [17, 18, 19, 20, 21]]
        #ethyl groups each including hydrogens
        e3 = [[0, 1, 7, 8, 9, 10, 11], [3, 4, 12, 13, 14, 15, 16],
              [5, 6, 17, 18, 19, 20, 21]]
        m1 = triethylamine.select(ethyl, hydrogen="exclude")
        m2 = triethylamine.select(ethyl, hydrogen="only")
        m3 = triethylamine.select(ethyl, hydrogen="include")
        self.assertEqual(e1, m1)
        self.assertEqual(e2, m2)
        self.assertEqual(e3, m3)

    def test_select_hydrogens(self):
        triethylamine = self.G.make_fragment("CCN(CC)CC")
        ethyl = "[#6][#6]"
        s = pybel.Smarts(ethyl)
        
        selected = s.findall(triethylamine.molecule)
        expectations = [[7, 8, 9, 10, 11],
                        [12, 13, 14, 15, 16],
                        [17, 18, 19, 20, 21]]

        for k, ethyl in enumerate(selected):
            #cinfony.pybel uses 0-based indexing while the lower level pybel
            #uses 1-based indexing, so subtract 1 from each index in match
            e = [s - 1 for s in ethyl]
            hydrogens = triethylamine.select_hydrogens(e)
            self.assertEqual(expectations[k], hydrogens)
            for h in hydrogens:
                self.assertEqual(1, triethylamine.atoms[h].atomicnum)

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
