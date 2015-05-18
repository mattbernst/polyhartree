#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_pdynamo
    ~~~~~~~~~~~~~~

    Test pDynamo specific functionality that is not handled elsewhere.
"""
import sys
import unittest
import geoprep
from tests.common_testcode import runSuite
from adapters import pdynamo

class PDTestCase(unittest.TestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = pdynamo.PDynamo()

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(PDTestCase, name = test_name)

    else:
        result = runSuite(PDTestCase)

    return result

if __name__ == "__main__":
    runTests()
