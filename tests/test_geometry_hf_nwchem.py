#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_geometry_hf_nwchem
    ~~~~~~~~~~~~~~

    Test NWChem implementations for Hartree-Fock geometry jobs.
"""
import sys
import geoprep
from adapters import nwchem
from tests import geometry_hf as gh
from tests.common_testcode import runSuite

class NWChemHFGeometryTestCase(gh.HFGeometryTestCase):
    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = nwchem.NWChem()

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(NWChemHFGeometryTestCase, name=test_name)

    else:
        result = runSuite(NWChemHFGeometryTestCase)

    return result

if __name__ == '__main__':
    runTests()
