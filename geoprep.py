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
        #iupac = webel.Molecule(molecule).write('iupac')
        iupac = "NOT A REAL IUPAC NAME"
        smiles = molecule.write('smiles')
        molecule.make3D()
        molecule.title = "{0} SMILES: {1}".format(iupac, smiles)
        enriched = self.enrich_molecule(molecule)
        return enriched

    
        
        
