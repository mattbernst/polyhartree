# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-

class MolecularCalculator(object):
    def create_geometry(self, molecule, options={}):
        raise NotImplementedError

    def extract_geometry(self, data, options={}):
        raise NotImplementedError

    def extract_energy(self, data, options={}):
        raise NotImplementedError

    def make_energy_job(self, molecule, method, options={}):
        raise NotImplementedError

    def make_opt_job(self, molecule, method, options={}):
        raise NotImplementedError

    def check_method(self, method):
        if method not in self.methods:
            raise ValueError("Unrecognized method {0}".format(repr(method)))

    def check_coordinates(self, coordinate_choice):
        if not coordinate_choice:
            return 'cartesian'

        elif coordinate_choice in ['zmatrix', 'cartesian']:
            return coordinate_choice

        else:
            raise ValueError("Unrecognized coordinates option {0}".format(repr(coordinate_choice)))
