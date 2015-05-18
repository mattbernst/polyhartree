#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_geometry_hf_psi4
    ~~~~~~~~~~~~~~

    Test Psi4 implementations for Hartree-Fock geometry jobs.
"""
import sys
import geoprep
from adapters import psi4
from tests.common_testcode import runSuite
from tests import geometry_hf as gh

class Psi4HFGeometryTestCase(gh.HFGeometryTestCase):
    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = psi4.Psi4()

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(Psi4HFGeometryTestCase, name=test_name)

    else:
        result = runSuite(Psi4HFGeometryTestCase)

    return result

if __name__ == '__main__':
    runTests()
