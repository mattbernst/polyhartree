#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_energy_hf_gamess_us
    ~~~~~~~~~~~~~~

    Test GAMESS-US implementations for Hartree-Fock energy jobs.
"""
import sys
import geoprep
from adapters import gamess_us
from tests.common_testcode import runSuite
from tests import energy_hf as eh

class GAMESSHFEnergyTestCase(eh.HFEnergyTestCase):
    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = gamess_us.GAMESSUS()

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(GAMESSHFEnergyTestCase, name = test_name)

    else:
        result = runSuite(GAMESSHFEnergyTestCase)

    return result

if __name__ == '__main__':
    runTests()
