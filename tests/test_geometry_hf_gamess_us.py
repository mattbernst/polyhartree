#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_geometry_hf_gamess_us
    ~~~~~~~~~~~~~~

    Test GAMESS-US implementations for Hartree-Fock geometry jobs.
"""

import sys
import geoprep
from adapters import gamess_us
from tests import geometry_hf as gh
from tests.common_testcode import runSuite

class GAMESSHFGeometryTestCase(gh.HFGeometryTestCase):
    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = gamess_us.GAMESSUS()

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(GAMESSHFGeometryTestCase, name=test_name)

    else:
        result = runSuite(GAMESSHFGeometryTestCase)

    return result

if __name__ == '__main__':
    runTests()
