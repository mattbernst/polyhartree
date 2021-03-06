# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_semiempirical_energies
    ~~~~~~~~~~~~~~

    Test energies using different implementations of common methods with the
    same geometries. This test was prompted by apparent problems with
    pdynamo energies for the methyl radical.
"""
import sys
import unittest
import geoprep
from tests.common_testcode import runSuite
from adapters import gamess_us
from adapters import mopac7
from adapters import pdynamo

class PDTestCase(unittest.TestCase):

    def setUp(self):
        self.G = geoprep.Geotool()
        self.implementations = {"pdynamo" : pdynamo.PDynamo(),
                                "gamess" : gamess_us.GAMESSUS(),
                                "mopac7" : mopac7.Mopac7()}

    def run_multi(self, geometry, method):
        """Run multiple energy calculations at the same geometry and method
        settings using different quantum chemistry back ends.

        :param geometry: system geometry for calculations
        :type geometry : cinfony.pybel.Molecule
        :param method: method name to use
        :type method : str
        :return: implementation and job results from each back end
        :rtype : list
        """
        
        results = []
        for implementation in sorted(self.implementations.keys()):
            k = self.implementations[implementation]
            job = k.make_energy_job(geometry, method)
            job.run()

            result = {"job" : job, "implementation" : implementation}
            results.append(result)

        return results

    def find_differences(self, results, metric):
        """Determine the differences in results. Use the first calculation
        in the list as reference value and see how much other results
        deviate from it.

        :param results: results of calculations with different back ends
        :type results : list
        :param metric: a measurable quantity to compare, e.g. "energy" or "hof"
        :type metric : str
        :return: calculation metric differences
        :rtype : list
        """

        def places_of_agreement(a, b):
            max_digits = 10
            for k in range(max_digits):
                try:
                    self.assertAlmostEqual(a, b, places=k)
                except AssertionError:
                    return k

            return max_digits
        
        transformed = []
        
        reference = getattr(results[0]["job"], metric)

        for result in results:
            current = getattr(result["job"], metric)
            delta = reference - current
            places = places_of_agreement(reference, current)

            changes = {"delta_{0}".format(metric) : delta,
                       "places_{0}".format(metric) : places}
            t = result.copy()
            t.update(changes)
            transformed.append(t)

        return transformed

    def _test_energy_differences(self, geometry, method, metric, min_places):
        """Compare a metric when running the same energy method with the same
        geometry across different chemistry back-ends.

        :param geometry: a trial geometry
        :type geometry : geoprep.System
        :param method: a method name, e.g. "semiempirical:am1"
        :type method : str
        :param metric: a measurable quantity to compare, e.g. "energy" or "hof"
        :type metric : str
        :param min_places: minimum number of decimal places agreement for metric
        :return: differences among packages, using the first run as reference
        :rtype : list
        """
        
        jobs = self.run_multi(geometry, method)
        differences = self.find_differences(jobs, metric)
        places_name = "places_{0}".format(metric)
        places = [x[places_name] for x in differences]

        self.assertTrue(min_places <= min(places))

        return differences

    def test_methylium(self):
        geometry = self.G.make_system("[CH3+]")
        self._test_energy_differences(geometry, "semiempirical:pm3",
                                      "heat_of_formation", 6)
        self._test_energy_differences(geometry, "semiempirical:am1",
                                      "heat_of_formation", 5)
        self._test_energy_differences(geometry, "semiempirical:mndo",
                                      "heat_of_formation", 5)

    def test_methane(self):
        geometry = self.G.make_system("C")
        self._test_energy_differences(geometry, "semiempirical:pm3",
                                      "heat_of_formation", 6)
        self._test_energy_differences(geometry, "semiempirical:am1",
                                      "heat_of_formation", 5)
        self._test_energy_differences(geometry, "semiempirical:mndo",
                                      "heat_of_formation", 5)

    def test_carbanide(self):
        geometry = self.G.make_system("[CH3-]")
        self._test_energy_differences(geometry, "semiempirical:pm3",
                                      "heat_of_formation", 6)
        self._test_energy_differences(geometry, "semiempirical:am1",
                                      "heat_of_formation", 5)
        self._test_energy_differences(geometry, "semiempirical:mndo",
                                      "heat_of_formation", 5)

    def test_methyl_radical(self):
        geometry = self.G.make_system("[CH3]")
        pm3 = self._test_energy_differences(geometry, "semiempirical:pm3",
                                            "heat_of_formation", 6)
        am1 = self._test_energy_differences(geometry, "semiempirical:am1",
                                            "heat_of_formation", 5)
        mndo = self._test_energy_differences(geometry, "semiempirical:mndo",
                                             "heat_of_formation", 5)

    def test_nitric_oxide(self):
        geometry = self.G.make_system("[N]=O")
        self._test_energy_differences(geometry, "semiempirical:pm3",
                                      "heat_of_formation", 5)
        self._test_energy_differences(geometry, "semiempirical:am1",
                                      "heat_of_formation", 5)
        self._test_energy_differences(geometry, "semiempirical:mndo",
                                      "heat_of_formation", 4)

    def test_ammonia(self):
        geometry = self.G.make_system("N")
        self._test_energy_differences(geometry, "semiempirical:pm3",
                                      "heat_of_formation", 4)
        self._test_energy_differences(geometry, "semiempirical:am1",
                                      "heat_of_formation", 5)
        self._test_energy_differences(geometry, "semiempirical:mndo",
                                      "heat_of_formation", 5)

    def test_water_dimer(self):
        #run a very crude water sort of water dimer calculation, since we
        #don't have geometry loading yet (no actual H-bond here!)
        w1 = self.G.make_fragment("O")
        w2 = self.G.make_fragment("O")
        w2.translate((0, 2.0, 0))
        s = geoprep.System([w1, w2])

        self._test_energy_differences(s, "semiempirical:pm3",
                                      "heat_of_formation", 4)
        self._test_energy_differences(s, "semiempirical:am1",
                                      "heat_of_formation", 4)
        self._test_energy_differences(s, "semiempirical:mndo",
                                      "heat_of_formation", 4)

    def test_translation_sensitivity(self):
        #Translating a system should not meaningfully affect its calculated
        #energy. But with Mopac7 there is an effect when translating linear
        #molecules like NO, CO, H2...

        M = mopac7.Mopac7()
        P = pdynamo.PDynamo()
        G = gamess_us.GAMESSUS()
        
        results = {}
        CO = self.G.make_system("[C-]#[O+]")
        formaldehyde = self.G.make_system("C=O")
        
        #translate a little on the x axis each time
        dx = 0.25
        tx = 0.0
        for j in range(4):
            r1 = (self.run_multi(CO, "semiempirical:pm3"), "CO")
            r2 = (self.run_multi(formaldehyde, "semiempirical:pm3"),
                  "formaldehyde")
            for r, name in (r1, r2):
                for entry in r:
                    key = entry["implementation"] + "-" + name
                    hof = entry["job"].heat_of_formation
                    try:
                        results[key].append(hof)
                    except KeyError:
                        results[key] = [hof]
            
            v = (dx, 0, 0)
            tx += dx
            CO.fragments[0].translate(v)
            formaldehyde.fragments[0].translate(v)

        #For linear carbon monoxide, mopac7 produces a different result with
        #each x-axis translation. Otherwise all calculations behave as
        #expected: translation does not change energy
        for key in results:
            num_unique = len(set(results[key]))
            if key == "mopac7-CO":
                self.assertEqual(len(results[key]), num_unique)
            else:
                self.assertEqual(1, num_unique)

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(PDTestCase, name = test_name)

    else:
        result = runSuite(PDTestCase)

    return result

if __name__ == "__main__":
    runTests()
