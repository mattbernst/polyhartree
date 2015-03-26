# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import os
import types
from cinfony import pybel, webel

class Fragment(object):
    def __init__(self, molecule):
        self.molecule = molecule

        #Try to pass through these attributes/methods of the underlying
        #molecule as attributes/methods of the fragment
        for attr in ['addh', 'atoms', 'calcdesc', 'calcfp', 'charge',
                     'conformers', 'data', 'dim', 'draw', 'energy',
                     'exactmass', 'formula', 'localopt', 'make3D', 'molwt',
                     'removeh', 'sssr', 'title', 'write']:
            try:
                p = getattr(molecule, attr)
                setattr(self, attr, p)
            except AttributeError:
                pass

    @property
    def nelec(self):
        """Return number of electrons in system

        @return: number of electrons
        @rtype : int
        """

        total = sum([a.atomicnum for a in self.atoms]) - self.charge
        return total

    @property
    def spin(self):
        return self.molecule.spin

    @spin.setter
    def spin(self, s):
        """Set spin by using the setter on the underlying hidden OBMol.

        @param s: spin to set
        @type s : int
        """
        
        self.molecule.OBMol.SetTotalSpinMultiplicity(s)

    def set_zero_to_origin(self):
        """Set coordinates of atom 0 to the origin coordinates: 0, 0, 0
        Also translate other atoms to match

        This may make Mopac7 a little less finicky.
        """

        base = self.molecule.atoms[0].coords
        neg = [a * -1 for a in base]
        self.translate(neg)

    def translate(self, vec):
        """Translate every atom in molecule by the coordinates in vec.

        It seems like it should be possible to call the underlying
        molecule.OBMol.Translate method directly, but unable to pass correct
        data type from python.

        @param molecule: molecule to be translated
        @type molecule : cinfony.pybel.Molecule
        @param vec: x, y, z coordinates
        @type vec : list
        """

        for k in range(len(self.molecule.atoms)):
            coords = self.molecule.atoms[k].coords
            translated = []
            for j in range(3):
                t = coords[j] + vec[j]
                translated.append(t)

            self.molecule.atoms[k].OBAtom.SetVector(*translated)


class Geotool(object):
    def make_fragment(self, representation, kind="smiles"):
        """Take a representation of a molecule and add a title with IUPAC
        and SMILES designations. Convert a linear or 2D molecular specification
        to a 3D form.

        TODO: Allow IUPAC, trivial-name input. Cache webel calls.

        @param representation: linear molecule encoding
        @type representation : str
        @param kind: smiles, inchi, etc. (default smiles)
        @type kind : str
        @return: 3D-form molecular fragment
        @rtype: Fragment
        """

        molecule = pybel.readstring(kind, representation)
        #iupac = webel.Molecule(molecule).write('iupac').strip()
        iupac = "NOT A REAL IUPAC NAME"
        smiles = molecule.write("smi").strip()
        molecule.make3D()
        fragment = Fragment(molecule)
        fragment.set_zero_to_origin()
        fragment.title = "{0} SMILES: {1}".format(iupac, smiles)
        return fragment

    def read_mol(self, file_name, fmt=None):
        """Load a molecular structure from a file. Guess at the file type
        from extension caller does not supply explicit fmt.

        TODO: reduce code duplication with make_mol

        @param file_name: file to open
        @type file_name : str
        @param fmt: optional OpenBabel format code e.g. "xyz"
        @type fmt : str
        @return: molecular system
        @rtype : cinfony.pybel.Molecule
        """

        if not fmt:
            try:
                fmt = file_name.split(".")[-1]
            except IndexError:
                msg = "No fmt given for {0} and unable to guess from file extensions".format(repr(file_name))

        molecule = pybel.readfile(fmt, file_name).next()
        iupac = "NOT A REAL IUPAC NAME"
        smiles = molecule.write("smi").strip()

        molecule.title = "{0} SMILES: {1}".format(iupac, smiles)
        enriched = self.enrich_molecule(molecule)
        return enriched
