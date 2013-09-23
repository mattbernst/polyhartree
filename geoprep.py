from cinfony import pybel, webel

class Geotool(object):
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
        iupac = webel.Molecule(molecule).write('iupac')
        smiles = molecule.write('smiles')
        molecule.make3D()
        molecule.title = "{0} SMILES: {1}".format(iupac, smiles)
        return molecule

    
        
        
