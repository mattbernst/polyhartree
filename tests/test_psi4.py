#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_psi4
    ~~~~~~~~~~~~~~

    Test Psi4 specific functionality that is not handled elsewhere.
"""

import sys
import unittest
import geoprep
from adapters import psi4
from tests.common_testcode import runSuite
from sharedutilities import Utility
from tests import reference_values

class Psi4TestCase(unittest.TestCase):
    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = psi4.Psi4()

    def test_make_control_block(self):
        #test control block formatting
        expected = "molecule {\n  0 1\n  H1 0.0 0.0 0.0\n  H2 0.9 0.0 0.0\n}"
        cartesians = ["molecule", "0 1", "H1 0.0 0.0 0.0", "H2 0.9 0.0 0.0"]
        block = self.C.make_control_block(cartesians)
        self.assertEqual(expected, block)

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
