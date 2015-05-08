#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_geometry_semiempirical_mopac7
    ~~~~~~~~~~~~~~

    Test MOPAC 7 semiempirical implementation for geometry.
"""

import sys
import geoprep
from adapters import mopac7
from tests.common_testcode import runSuite
from tests import geometry_semiempirical as gs
from tests import reference_values

class MOPACGeometryTestCase(gs.SemiempiricalGeometryTestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = mopac7.Mopac7()

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(MOPACGeometryTestCase, name = test_name)

    else:
        result = runSuite(MOPACGeometryTestCase)

    return result

if __name__ == "__main__":
    runTests()
