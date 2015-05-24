#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_gamess
    ~~~~~~~~~~~~~~

    Test GAMESS-US specific functionality that is not handled elsewhere.
"""
import sys
import unittest
import geoprep
from adapters import gamess_us
from tests import adapter
from tests import reference_values
from tests.common_testcode import runSuite

class GAMESSTestCase(adapter.AdapterTestCase):
    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = gamess_us.GAMESSUS()

    def test_create_geometry_explicit_symmetry(self):
        #test explicit symmetry setting (creation) for GAMESS-US
        #N.B.: actual specified geometry is wrong for TD -- job will not work!
        methane = self.G.make_fragment("C")
        methane.set_basis_name("3-21G")
        symmetry_group = "TD"
        options = {"symmetry" : symmetry_group}
        job = self.C.make_energy_job(methane, "hf:rhf", options=options)
        self.assertTrue(symmetry_group in job.deck)

    def test_bad_input_error(self):
        job = self.get_job("C", "semiempirical:pm3")

        #introduce an error in the input deck: misspell SCFTYP as SCFTYPE
        job.deck = job.deck.replace("SCFTYP", "SCFTYPE")

        self.assertEqual("begin", job.runstate)
        job.run()

        self.assertEqual("error", job.runstate)

    def test_reformat_long_line(self):
        #verify splitting of long directives so as to comfortably fit within
        #the 80 character per line limit of GAMESS-US
        maxlen = 50
        line = " $BASIS basnam(1)=Cl0, As0, Cl0, C0, C0, Cl0, H0, Cl1 $END"
        pieces, begin, end = self.C.reformat_long_line(line, " $BASIS", "$END",
                                                       maxlen=maxlen)
        for piece in pieces:
            self.assertTrue(len(piece) <= maxlen)

    def test_prepare_basis_data(self):
        #test generation of inline basis set data with hydrogen peroxide
        #should generate one basis set assignment for oxygen and two for
        #hydrogen, since one hydrogen gets assigned a different basis
        btag = "basis_tag"
        expected_labels = ["$H0", "$H1", "$O0"]
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

    def test_extract_geometry_from_log(self):
        #read/verify geometry from a specific stored RHF water optimization
        self.geometry_from_log("O", "hf:rhf",
                               "tests/data/logs/rhf_water_geo_min_gamess_us.log",
                               6)

    def test_extract_energy_from_log(self):
        #read/verify geometry from a specific stored single point energy log
        self.energy_from_log("CO", "hf:rhf",
                             "tests/data/logs/rhf_methanol_energy_gamess_us.log",
                             reference_values.methanol_rhf_321g, None)

    def test_extract_hof_from_log(self):
        #read/verify geometry from a specific stored single point energy log
        self.energy_from_log("C", "semiempirical:pm3",
                             "tests/data/logs/pm3_methane_energy_gamess_us.log",
                             None, reference_values.methane_pm3_hof)



def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(GAMESSTestCase, name = test_name)

    else:
        result = runSuite(GAMESSTestCase)

    return result

if __name__ == '__main__':
    runTests()
