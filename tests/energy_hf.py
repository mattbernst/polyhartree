#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    energy_hf
    ~~~~~~~~~~~~~~

    Test common Hartree-Fock implementations for energy jobs.
"""

import sys
import unittest
import geoprep
from tests import reference_values
from tests import energy

class HFEnergyTestCase(energy.EnergyTestCase):
    def test_energy_scf_methane(self, places=6):
        #very basic minimal basis set test for methane
        methane = self.G.make_fragment("C")
        methane.set_basis_name("3-21G")
        job = self.C.make_energy_job(methane, "hf:rhf")
        job.run()
        self.assertNearMatch(reference_values.methane_rhf_321g,
                             job.energy, places=places)

    def test_energy_rohf_uhf_scf_methane(self, places=6):
        #compare UHF and ROHF across methyl radicals for HF energy
        
        smiles = {"methyl_radical" : "[CH3]"}
        jobs = {}

        for key in smiles:
            s = smiles[key]
            fragment = self.G.make_fragment(s)
            fragment.set_basis_name("cc-pVDZ")

            for reference in ["uhf", "rohf"]:
                j = self.C.make_energy_job(fragment,
                                           "hf:{0}".format(reference))
                self.assertTrue(reference in j.deck.lower())
                j.run()
                name = "{0}-{1}".format(key, reference)
                jobs[name] = j

        energies = dict([(j, jobs[j].energy) for j in jobs])

        self.assertNearMatch(reference_values.methyl_rohf_ccpvdz,
                             energies["methyl_radical-rohf"], places=places)
        self.assertNearMatch(reference_values.methyl_uhf_ccpvdz,
                             energies["methyl_radical-uhf"], places=places)

    def test_energy_methanol(self, places=6):
        #regression test using geometry read from external file
        methanol = self.G.read_fragment("tests/data/methanol1.xyz")
        methanol.set_basis_name("3-21G")
        job = self.C.make_energy_job(methanol, "hf:rhf")
        job.run()
        self.assertNearMatch(reference_values.methanol_rhf_321g,
                             job.energy, places=places)
