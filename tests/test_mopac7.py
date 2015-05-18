#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_mopac
    ~~~~~~~~~~~~~~

    Test MOPAC 7 specific functionality that is not handled elsewhere.
"""
import sys
import unittest
import geoprep
from tests import adapter
from tests.common_testcode import runSuite
from adapters import mopac7

class MOPACTestCase(adapter.AdapterTestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = mopac7.Mopac7()

    def test_bad_input_error(self):
        #This fails because we introduce a typo into the deck.
        job = self.get_job("C", "semiempirical:pm3")
        job.deck = job.deck.replace("SINGLET", "SNIGLET")
        job.run()
        self.assertEqual("error", job.runstate)

    def test_element_error(self):
        #Try to create a system containing an element unparameterized for
        #MINDO/3. Should raise an error.
        #N.B.: The generated input file would actually run in practice,
        #rather contrary to the Mopac 7 manual's parameterization information,
        #but it seems best to err on the side of caution.
        stibine = self.G.make_system("[SbH3]")
        self.assertRaises(ValueError, self.C.make_energy_job, stibine,
                          "semiempirical:mindo/3")

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(MOPACTestCase, name = test_name)

    else:
        result = runSuite(MOPACTestCase)

    return result

if __name__ == "__main__":
    runTests()
