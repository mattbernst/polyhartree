#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_geometry_semiempirical_gamess_us
    ~~~~~~~~~~~~~~

    Test GAMESS-US semiempirical implementation for geometry.
"""

import sys
import geoprep
from adapters import gamess_us
from tests.common_testcode import runSuite
from tests import geometry_semiempirical as gs
from tests import reference_values

class GAMESSGeometryTestCase(gs.SemiempiricalGeometryTestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = gamess_us.GAMESSUS()

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(GAMESSGeometryTestCase, name = test_name)

    else:
        result = runSuite(GAMESSGeometryTestCase)

    return result

if __name__ == "__main__":
    runTests()
