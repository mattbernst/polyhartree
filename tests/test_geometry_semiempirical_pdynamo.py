#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_geometry_semiempirical_pdynamo
    ~~~~~~~~~~~~~~

    Test pDynamo implementations for semiempirical geometry.
"""
import sys
import geoprep
from adapters import pdynamo
from tests.common_testcode import runSuite
from tests import geometry_semiempirical as gs

class PDGeometryTestCase(gs.SemiempiricalGeometryTestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = pdynamo.PDynamo()

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(PDGeometryTestCase, name = test_name)

    else:
        result = runSuite(PDGeometryTestCase)

    return result

if __name__ == "__main__":
    runTests()
