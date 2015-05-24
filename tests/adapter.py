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

    def read_job_log(self, smiles, theory, log_file):
        """Get contents from a previously run job log and a matching
        :param smiles: a SMILES linear representation of a molecule
        :type smiles : str
        :param theory: theory to use, e.g. "hf:rhf" or "semiempirical:pm3"
        :type theory : str
        :param basis: optional basis set name, for calculations using basis sets
        :type basis : None or str
        :return: a job appropriate to the log file, log file contents
        :rtype : dict
        """

        basis = None
        if "semiempirical" not in theory:
            basis = "3-21G"

        job = self.get_job(smiles, theory, basis=basis)
        with open(log_file) as infile:
            data = infile.read()

        r = {"data" : data, "job" : job}
        return r

    def test_extract_geometry_from_log(self, smiles, theory, log_file, n_snapshots):
        #read/verify geometry from a specific stored RHF water optimization
        components = self.read_job_log(smiles, theory, log_file)
        job, data = components["job"], components["data"]

        self.assertEqual(0, len(job.geometry_history))
        job.extract_geometry(data)
        self.assertEqual(n_snapshots, len(job.geometry_history))
