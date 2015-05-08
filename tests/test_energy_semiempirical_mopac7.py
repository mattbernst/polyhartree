#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_energy_semiempirical_mopac7
    ~~~~~~~~~~~~~~

    Test MOPAC 7 semiempirical implementation for energy.
"""

import sys
import geoprep
from adapters import mopac7
from tests.common_testcode import runSuite
from tests import energy_semiempirical as es
from tests import reference_values

class MOPACEnergyTestCase(es.SemiempiricalEnergyTestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = mopac7.Mopac7()

    def test_energy_mindo3_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:mindo/3")
        job.run()
        self.assertNearMatch(reference_values.methane_mindo3_hof,
                             job.heat_of_formation, places=5)

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(MOPACEnergyTestCase, name = test_name)

    else:
        result = runSuite(MOPACEnergyTestCase)

    return result

if __name__ == "__main__":
    runTests()
