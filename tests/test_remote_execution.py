#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_remote_execution
    ~~~~~~~~~~~~~~

    Test job execution over ansible, using ssh transport to 127.0.0.1. Use
    NWChem and MOPAC 7 as back ends for test purposes because they are the
    best supported via Debian packages.
"""

import sys
import unittest
import geoprep
from adapters import nwchem, mopac7
from sharedutilities import Utility

class RETestCase(unittest.TestCase):
    def setUp(self):
        self.G = geoprep.Geotool()
        self.C = nwchem.NWChem()
        self.C2 = mopac7.Mopac7()

    def tearDown(self):
        pass

    def test_nwchem_energy_scf_methane(self):
        #very basic minimal basis set test for methane
        expected_energy = -39.976642
        methane = self.G.make_fragment("C")
        methane.set_basis_name("3-21G")

        #local execution
        job = self.C.make_energy_job(methane, "hf:rhf")
        job.run()
        self.assertAlmostEqual(expected_energy, job.energy, places=5)

        #over ansible
        job2 = self.C.make_energy_job(methane, "hf:rhf")
        job2.run(host="127.0.0.1")
        self.assertAlmostEqual(expected_energy, job2.energy, places=5)

    def test_nwchem_bad_input_error(self):
        #job will terminate abnormally
        methane = self.G.make_fragment("C")
        methane.set_basis_name("3-21G")
        job = self.C.make_energy_job(methane, "hf:rhf")

        #introduce an error in the input deck: misspell rhf as thf
        job.deck = job.deck.replace("RHF", "THF")

        #local execution
        self.assertEqual("begin", job.runstate)
        job.run()
        self.assertEqual("error", job.runstate)

        #over ansible
        job2 = self.C.make_energy_job(methane, "hf:rhf")
        job2.deck = job.deck.replace("RHF", "THF")
        self.assertEqual("begin", job2.runstate)
        job2.run(host="127.0.0.1")
        self.assertEqual("error", job2.runstate)

    def test_mopac_energy_pm3_methylium(self):
        expected_energy = -5.641481
        expected_hof = 0.408868
        
        methylium = self.G.make_system("[CH3+]")

        #local execution
        job = self.C2.make_energy_job(methylium, "semiempirical:pm3")
        job.run()
        self.assertAlmostEqual(expected_energy, job.energy, places=5)
        self.assertAlmostEqual(expected_hof, job.heat_of_formation, places=5)

        #over ansible
        job2 = self.C2.make_energy_job(methylium, "semiempirical:pm3")
        job2.run(host="127.0.0.1")
        self.assertAlmostEqual(expected_energy, job2.energy, places=5)
        self.assertAlmostEqual(expected_hof, job2.heat_of_formation, places=5)

    def test_ansible_bad_host(self):
        #try to use ansible on a badly configured host and confirm error state
        methylium = self.G.make_system("[CH3+]")
        job = self.C2.make_energy_job(methylium, "semiempirical:pm3")
        result = job.ansible_run("shell", "ls", "notarealhost")
        self.assertEqual(1, len(job.messages))
        err = job.messages[0]
        self.assertTrue("unknown error" in err or "sshpass" in err)

    def test_ansible_missing_host(self):
        #try to use ansible on a non-configured host and get an exception
        methylium = self.G.make_system("[CH3+]")
        job = self.C2.make_energy_job(methylium, "semiempirical:pm3")
        self.assertRaises(KeyError, job.ansible_run, "shell", "ls",
                          "nosuchhost")

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
        result = runSuite(RETestCase, name = test_name)

    else:
        result = runSuite(RETestCase)

    return result

if __name__ == '__main__':
    runTests()
