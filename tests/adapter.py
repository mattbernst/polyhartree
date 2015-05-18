#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    adapter
    ~~~~~~~~~~~~~~

    Common code for adapter tests.
"""
import unittest

class AdapterTestCase(unittest.TestCase):
    def get_job(self, smiles, theory, basis=None):
        """Get a job we can use to test parts of the job code.

        :param smiles: a SMILES linear representation of a molecule
        :type smiles : str
        :param theory: theory to use, e.g. "hf:rhf" or "semiempirical:pm3"
        :type theory : str
        :param basis: optional basis set name, for calculations using basis sets
        :type basis : None or str
        :return: an energy job
        """

        m = self.G.make_fragment(smiles)
        if basis:
            m.set_basis_name(basis)
        job = self.C.make_energy_job(m, theory)

        return job
