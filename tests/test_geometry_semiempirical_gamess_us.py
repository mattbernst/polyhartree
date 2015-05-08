#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_geometry_semiempirical_gamess_us
    ~~~~~~~~~~~~~~

    Test GAMESS-US semiempirical implementation for geometry.
"""

import sys
import unittest
import geoprep
from adapters import gamess_us
from tests import geometry_semiempirical as gs
from tests import reference_values

class GAMESSGeometryTestCase(gs.SemiempiricalGeometryTestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = gamess_us.GAMESSUS()

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
        result = runSuite(GAMESSGeometryTestCase, name = test_name)

    else:
        result = runSuite(GAMESSGeometryTestCase)

    return result

if __name__ == "__main__":
    runTests()
