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
from tests import adapter
from tests.common_testcode import runSuite
from adapters import pdynamo

class PDTestCase(adapter.AdapterTestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = pdynamo.PDynamo()

    def test_extract_geometry_from_log(self):
        #read/verify geometry from a specific stored PM3 water optimization
        job = self.get_job("O", "semiempirical:pm3")
        with open("tests/data/logs/pm3_water_geo_min_pdynamo.log") as infile:
            data = infile.read()

        self.assertEqual(0, len(job.geometry_history))
        job.extract_geometry(data)
        self.assertEqual(9, len(job.geometry_history))

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
