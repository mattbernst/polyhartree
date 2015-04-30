#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_psi4
    ~~~~~~~~~~~~~~

    Test Psi4 implementations for energy, minimum geometry,
    transition state geometry, and frequencies.
"""

import sys
import unittest
import geoprep
from adapters import psi4
from sharedutilities import Utility

class Psi4TestCase(unittest.TestCase):
    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = psi4.Psi4()

    def tearDown(self):
        pass

    def test_make_control_block(self):
        #test control block formatting
        expected = "molecule {\n  0 1\n  H1 0.0 0.0 0.0\n  H2 0.9 0.0 0.0\n}"
        cartesians = ["molecule", "0 1", "H1 0.0 0.0 0.0", "H2 0.9 0.0 0.0"]
        block = self.C.make_control_block(cartesians)
        self.assertEqual(expected, block)

    def test_energy_scf_methane(self):
        #very basic minimal basis set test for methane
        methane = self.G.make_fragment("C")
        methane.set_basis_name("3-21G")
        job = self.C.make_energy_job(methane, "hf:rhf")
        job.run()
        self.assertAlmostEqual(-39.976642, job.energy, places=5)

    def test_energy_rohf_uhf_scf_methane(self):
        #compare UHF and ROHF across methyl radicals for HF energy
        
        smiles = {"methyl_radical" : "[CH3]"}
        jobs = {}

        for key in smiles:
            s = smiles[key]
            fragment = self.G.make_fragment(s)
            fragment.set_basis_name("cc-pVDZ")

            for reference in ["uhf", "rohf"]:
                j = self.C.make_energy_job(fragment,
                                           "hf:{0}".format(reference))
                self.assertTrue(reference in j.deck.lower())
                j.run()
                name = "{0}-{1}".format(key, reference)
                jobs[name] = j

        energies = dict([(j, jobs[j].energy) for j in jobs])

        #These match NWChem/GAMESS to only 4 decimal places due to default
        #use of density fitting. Direct SCF matches to 8 places.
        self.assertAlmostEqual(-39.5596229, energies["methyl_radical-rohf"],
                               places=4)
        self.assertAlmostEqual(-39.5638132, energies["methyl_radical-uhf"],
                               places=4)

    def test_bad_input_error(self):
        methane = self.G.make_fragment("C")
        methane.set_basis_name("3-21G")
        job = self.C.make_energy_job(methane, "hf:rhf")

        #introduce an error in the input deck: misspell scf as sfc
        job.deck = job.deck.replace("scf", "sfc")

        self.assertEqual("begin", job.runstate)
        job.run()

        self.assertEqual("error", job.runstate)

    def test_prepare_basis_data(self):
        #test generation of inline basis set data with hydrogen peroxide
        #should generate one basis set assignment for oxygen and two for
        #hydrogen, since one hydrogen gets assigned a different basis
        btag = "basis_tag"
        expected_labels = ["H0", "H1", "O0"]
        expected_tags = ["O0", "O0", "H0", "H1"]

        peroxide = self.G.make_fragment("OO")
        peroxide.set_basis_name("3-21G")
        hydrogens = peroxide.select("[O]", hydrogen="only")
        peroxide.set_basis_name("6-31G", hydrogens[:1])
        
        s = geoprep.System([peroxide])
        b = self.C.prepare_basis_data(s, options={"basis_tag_name" : btag})

        for label in expected_labels:
            self.assertTrue(label in b["basis_data"])

        basis_tags = s.atom_properties(btag)
        self.assertEqual(expected_tags, basis_tags)

    def test_make_tagged_cartesian_geometry(self):
        #test generation of cartesian coordinates matched to atom tags
        btag = "basis_tag"
        options = {"basis_tag_name" : btag, "property_name" : btag}
        expected_tags = ["O0", "O0", "H0", "H1"]

        peroxide = self.G.make_fragment("OO")
        peroxide.set_basis_name("3-21G")
        hydrogens = peroxide.select("[O]", hydrogen="only")
        peroxide.set_basis_name("6-31G", hydrogens[:1])
        
        s = geoprep.System([peroxide])
        self.C.prepare_basis_data(s, options=options)
        coordinates = self.C.make_tagged_cartesian_geometry(s, options=options)

        u = Utility()
        for j, c in enumerate(coordinates):
            self.assertTrue(c.startswith(expected_tags[j]))
            numbers = u.numericize(c, numeric_only=True)
            self.assertEqual(3, len(numbers))

    def test_create_geometry_basic(self):
        #test generation of complete geometry block
        btag = "basis_tag"
        expected_tags = ["O0", "H0"]
        options = {"basis_tag_name" : btag, "property_name" : btag}

        peroxide = self.G.make_fragment("OO")
        peroxide.set_basis_name("3-21G")
        
        s = geoprep.System([peroxide])
        self.C.prepare_basis_data(s, options=options)
        geometry = self.C.create_geometry(s, options=options)

        for tag in expected_tags:
            self.assertTrue(tag in geometry)

    def test_create_geometry_c1_symmetry(self):
        #test generation of geometry block with explicit C1 symmetry
        btag = "basis_tag"
        options = {"basis_tag_name" : btag, "property_name" : btag,
                   "symmetry" : "C1"}

        peroxide = self.G.make_fragment("OO")
        peroxide.set_basis_name("3-21G")
        
        s = geoprep.System([peroxide])
        self.C.prepare_basis_data(s, options=options)
        geometry = self.C.create_geometry(s, options=options)

        self.assertTrue("symmetry C1" in geometry)

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
        result = runSuite(Psi4TestCase, name = test_name)

    else:
        result = runSuite(Psi4TestCase)

    return result

if __name__ == '__main__':
    runTests()
