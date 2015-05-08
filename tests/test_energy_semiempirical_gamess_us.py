#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_energy_semiempirical_gamess_us
    ~~~~~~~~~~~~~~

    Test GAMESS-US implementations for semiempirical energy jobs.
"""

import sys
import unittest
import geoprep
from adapters import gamess_us
from tests import energy_semiempirical as es
from tests import reference_values

class GAMESSTestCase(es.SemiempiricalEnergyTestCase):
    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = gamess_us.GAMESSUS()

    def tearDown(self):
        pass

    def test_energy_rohf_uhf_pm3(self):
        #compare UHF and ROHF across different radicals for heat of formation
        
        smiles = {"methyl_radical" : "[CH3]",
                  "phenyl_radical" : "c1cc[c]cc1"}
        jobs = {}

        for key in smiles:
            s = smiles[key]
            mol = self.G.make_system(s)
            for reference in ["uhf", "rohf"]:
                j = self.C.make_energy_job(mol, "semiempirical:pm3",
                                           options={"reference" : reference})
                j.run()
                name = "{0}-{1}".format(key, reference)
                jobs[name] = j

        hofs = dict([(j, jobs[j].heat_of_formation) for j in jobs])

        self.assertNearMatch(reference_values.phenyl_rohf_pm3_hof,
                             hofs["phenyl_radical-rohf"], places=5)
        self.assertNearMatch(reference_values.phenyl_uhf_pm3_hof,
                             hofs["phenyl_radical-uhf"], places=5)
        self.assertNearMatch(reference_values.methyl_rohf_pm3_hof,
                             hofs["methyl_radical-rohf"], places=5)
        self.assertNearMatch(reference_values.methyl_uhf_pm3_hof,
                             hofs["methyl_radical-uhf"], places=5)


    def test_energy_rm1_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:rm1")
        job.run()
        self.assertNearMatch(reference_values.methane_rm1_hof,
                             job.heat_of_formation, places=5)

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
        result = runSuite(GAMESSTestCase, name = test_name)

    else:
        result = runSuite(GAMESSTestCase)

    return result

if __name__ == '__main__':
    runTests()
