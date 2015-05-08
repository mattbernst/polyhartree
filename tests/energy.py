#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    energy
    ~~~~~~~~~~~~~~

    Common energy job testing.
"""

import sys
import unittest
import geoprep
from tests import reference_values

class EnergyTestCase(unittest.TestCase):

    def assertNearMatch(self, reference_value, trial_value, places=6):
        """Check that the reference value and the trial value match to at
        least places.

        @param reference_value: value treated as reference for comparison
        @type reference_value : float
        @param trial_value: value tested against reference
        @type trial_value : float
        @param places: number of decimal places to match reference and trial
        @type places : int
        """

        self.assertAlmostEqual(reference_value, trial_value, places)
