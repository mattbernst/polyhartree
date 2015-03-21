import os
import types
from cinfony import pybel, webel

class Geotool(object):
    def enrich_molecule(self, molecule):
        """Monkey patch a pybel molecule to add new methods.

        @param molecule: molecule to gain new methods
        @type molecule : cinfony.pybel.Molecule
        @return: enriched molecule
        @rtype : cinfony.pybel.Molecule
        """

        def nelec(self):
            """Return number of electrons in system

            @return: number of electrons
            @rtype : int
            """

            total = sum([a.atomicnum for a in self.atoms]) + self.charge
            return total

        molecule.nelec = types.MethodType(nelec, molecule)
        return molecule

    def set_zero_to_origin(self, molecule):
        """Set coordinates of atom 0 to the origin coordinates: 0, 0, 0
        Also translate other atoms to match

        This may make Mopac7 a little less finicky.

        @param molecule: a molecule with 3D geometry
        @type molecule: cinfony.pybel.Molecule
        """

        base = molecule.atoms[0].coords
        neg = [a * -1 for a in base]
        self.translate(molecule, neg)

    def translate(self, molecule, vec):
        """Translate every atom in molecule by the coordinates in vec.

        It seems like it should be possible to call the underlying
        molecule.OBMol.Translate method directly, but unable to pass correct
        data type from python.

        @param molecule: molecule to be translated
        @type molecule : cinfony.pybel.Molecule
        @param vec: x, y, z coordinates
        @type vec : list
        """

        for k in range(len(molecule.atoms)):
            coords = molecule.atoms[k].coords
            translated = []
            for j in range(3):
                t = coords[j] + vec[j]
                translated.append(t)

            molecule.atoms[k].OBAtom.SetVector(*translated)
            
    def make_mol(self, representation, kind="smiles"):
        """Take a representation of a molecule and add a title with IUPAC
        and SMILES designations. Convert a linear or 2D molecular specification
        to a 3D form.

        TODO: Allow IUPAC, trivial-name input. Cache webel calls.

        @param representation: linear molecule encoding
        @type representation : str
        @param kind: smiles, inchi, etc. (default smiles)
        @type kind : str
        @return: 3D-form molecule
        @rtype: cinfony.pybel.Molecule
        """

        molecule = pybel.readstring(kind, representation)
        #iupac = webel.Molecule(molecule).write('iupac').strip()
        iupac = "NOT A REAL IUPAC NAME"
        smiles = molecule.write("smi").strip()
        molecule.make3D()
        self.set_zero_to_origin(molecule)
        molecule.title = "{0} SMILES: {1}".format(iupac, smiles)
        enriched = self.enrich_molecule(molecule)
        return enriched

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
        

    
        
        
