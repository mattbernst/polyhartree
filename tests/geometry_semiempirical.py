#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    geometry_semiempirical
    ~~~~~~~~~~~~~~

    Test common semiempirical implementations for geometry jobs.
"""

import sys
import unittest
import geoprep

class SemiempiricalGeometryTestCase(unittest.TestCase):
    def test_extract_geometry_unchanged(self):
        #extract geometry from an energy job, which should be basically the
        #same as the input
        methane = self.G.make_system("C")
        job = self.C.make_energy_job(methane, "semiempirical:pm3")
        job.run()

        reread = self.G.geolist_to_fragment(job.geometry)
        rmsd = self.G.align(methane.fragments[0], reread)["rmsd"]
        self.assertTrue(rmsd < 10**-4)

    def test_opt_pm3_water(self):
        #minimize geometry and verify that final structure has lower energy
        water = self.G.make_system("O")
        job = self.C.make_opt_job(water, "semiempirical:pm3")
        job.run()

        self.assertTrue(len(job.geometry_history) > 1)
        
        first = self.G.geolist_to_fragment(job.geometry_history[0])
        job1 = self.C.make_energy_job(first, "semiempirical:pm3")
        job1.run()

        self.assertTrue(job.energy < job1.energy)
