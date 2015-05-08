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
        self.assertTrue(rmsd < 10**-5)
