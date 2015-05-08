#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    energy_semiempirical
    ~~~~~~~~~~~~~~

    Test common semiempirical implementations for energy jobs.
"""

import sys
import unittest
import geoprep
from tests import reference_values

class SemiempiricalEnergyTestCase(unittest.TestCase):

    def setUp(self):
        self.G = geoprep.Geotool()

    def test_energy_pm3_methylium(self):
        methylium = self.G.make_system("[CH3+]")
        job = self.C.make_energy_job(methylium, "semiempirical:pm3")
        job.run()
        self.assertAlmostEqual(reference_values.methylium_pm3_hof,
                               job.heat_of_formation, places=5)

    def test_energy_pm3_carbanide(self):
        carbanide = self.G.make_system("[CH3-]")
        job = self.C.make_energy_job(carbanide, "semiempirical:pm3")
        job.run()
        self.assertAlmostEqual(reference_values.carbanide_pm3_hof,
                               job.heat_of_formation, places=5)

    def test_energy_pm3_methyl_radical(self):
        mradical = self.G.make_system("[CH3]")
        job = self.C.make_energy_job(mradical, "semiempirical:pm3")
        self.assertEqual("Forcing UHF for multiplicity 2", self.C.messages[0])
        job.run()
        self.assertAlmostEqual(reference_values.methyl_uhf_pm3_hof,
                               job.heat_of_formation, places=5)

    def test_energy_pm3_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pm3")
        job.run()
        self.assertAlmostEqual(reference_values.methane_pm3_hof,
                               job.heat_of_formation, places=5)

    def test_energy_mndo_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:mndo")
        job.run()
        self.assertAlmostEqual(reference_values.methane_mndo_hof,
                               job.heat_of_formation, places=5)

    def test_energy_am1_methane(self):
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:am1")
        job.run()
        self.assertAlmostEqual(reference_values.methane_am1_hof,
                               job.heat_of_formation, places=5)
