# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import os
import string
from cinfony import pybel, webel

class System(object):
    def __init__(self, fragment_or_fragments, spin=None, title=None):
        """Create a System containing one or more fragments. If there is just
        one fragment it will be promoted to a one-item list; sequences will
        become lists.

        @param fragment_or_fragments: a fragment or sequence of fragments
        @type fragment_or_fragments : Fragment | tuple | list
        @param spin: optional spin multiplicity for the entire system
        @type spin : int
        @param title: optional title for the entire system
        @type title : str
        """
        
        try:
            self.fragments = [f for f in fragment_or_fragments]
        except TypeError:
            self.fragments = [fragment_or_fragments]

        self.explicit_spin = spin
        self.explicit_title = title

    @property
    def charge(self):
        """Return sum of charges across all fragments to give system charge

        @return: total charge
        @rtype : int
        """

        total = sum([f.charge for f in self.fragments])
        return total

    @property
    def nelec(self):
        """Return number of electrons in system

        @return: number of electrons
        @rtype : int
        """

        total = sum([f.nelec for f in self.fragments])
        return total

    @property
    def title(self):
        """Return the explicitly set title, if present, or the title of the
        first fragment otherwise.

        @return: title
        @rtype : str
        """

        t = self.explicit_title
        if t is None:
            t = self.fragments[0].title

        return t

    @title.setter
    def title(self, title):
        """Set the explicit title.

        @param title: new title to set
        @type title : str
        """

        self.explicit_title = title

    @property
    def spin(self):
        """Return the explicitly set spin multiplicity, if present, or the
        spin multiplicity of the first fragment otherwise. This is probably
        an ok guess if you have all closed shell fragments or just one open
        shell fragment, appearing first in the list.

        N.B.: Assigning/computing with the "right" multiplicity is non-trivial
        for transition states, open shell molecule clusters, diradicals,
        transition metals, and probably more.

        @return: spin multiplicity for the system
        @rtype : int
        """
        
        if self.explicit_spin is not None:
            s = self.explicit_spin
        else:
            s = self.fragments[0].spin

        return s

    @spin.setter
    def spin(self, s):
        """Set explicit spin multiplicity

        @param s: spin to set
        @type s : int
        """
        
        self.explicit_spin = s

    @property
    def atoms(self):
        """Get all atoms in a system as a concatenated list of atoms in each
        fragment in the system.

        @return: all system atoms
        @rtype : list
        """

        atoms = []
        for f in self.fragments:
            atoms += f.atoms

        return atoms

    @property
    def elements(self):
        """Get all elements incorporated in system.

        @param system: System
        @type system : geoprep.System
        @return: unique list of system elements, in order of atomic number
        @rtype : list
        """
        
        elements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
                    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
                    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
                    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
                    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
                    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
                    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
                    "Md", "No", "Lr"]

        included_elements = set()

        for atom in self.atoms:
            symbol = elements[atom.atomicnum - 1]
            included_elements.add(symbol)

        indexed = [(elements.index(e), e) for e in included_elements]
        indexed.sort()
        finals = [e[1] for e in indexed]
        return finals

    def write(self, fmt):
        """Write all fragments into a single molecular system, by way of an
        intermediate fused xyz. There might be a better way to do this,
        but it's not obvious.

        There may be strange results trying to write systems containing
        multiple fragments as a linear format like SMILES or InChI

        @param fmt: final format to write, e.g. "mopin" or "pdb"
        @type fmt : str
        """

        natoms = len(self.atoms)
        xyzs = []
        for f in self.fragments:
            #clip header and get right to the atom geometry
            xyz = f.write("xyz").strip().split("\n\n")[1]
            xyzs.append(xyz)

        geoblock = "\n".join(xyzs)
        fused = "{0}\n\n{1}\n".format(natoms, geoblock)
        reread = pybel.readstring("xyz", fused)
        reread.title = self.title

        written = reread.write(fmt)

        return written


class Fragment(object):
    def __init__(self, molecule):
        """Create a Fragment which is a wrapper around one of the cinfony
        Molecule types (currently pybel only tested) with extra convenience
        methods.

        @param molecule: underlying cinfony molecule
        @type molecule : cinfony.*.Molecule
        """
        
        #keep a "universal SMILES" representation of the fragment handy
        m = pybel.Molecule(molecule)
        self.smiles = m.write("smi", opt={"U" : True})
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
        """Return number of electrons in fragment

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

    def select(self, smarts, hydrogen="include"):
        """Select atoms matching a SMARTS pattern with different treatments
        for hydrogen atoms.

        http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

        Hydrogen treatment:
        "include": hydrogen atoms attached to heavies become part of selection
        "exclude": return heavy atoms only
        "only": hydrogen atoms attached to heavies are returned without heavies

        @param smarts: a SMARTS pattern
        @type smarts : str
        @param hydrogen: "include", "exclude", or "only" (see above)
        @type hydrogen : str
        @return: indexes of atoms matching the selection
        @rtype : list
        """

        finder = pybel.Smarts(smarts)
        h = finder.findall(self.molecule)
        heavies = []
        results = []

        #shift indexes of initial heavy atom match to compensate for
        #0-based indexing in cinfony.pybel
        for group in h:
            g = [i - 1 for i in group]
            heavies.append(g)

        if hydrogen == "exclude":
            results = heavies

        elif hydrogen == "only":
            for group in heavies:
                hydrogens = self.select_hydrogens(group)
                results.append(hydrogens)

        elif hydrogen == "include":
            for group in heavies:
                hydrogens = self.select_hydrogens(group)
                merged = group + hydrogens
                results.append(merged)

        else:
            raise ValueError("Unrecognized option for hydrogen selection")
            
        return results

    def select_hydrogens(self, atoms):
        """Get indexes of any hydrogen atoms attached to input atom indexes.

        N.B.: Internal pybel indexes use 1-based indexing and cinfony.pybel
        uses 0-based indexing, so indexes of target atoms need to be reduced
        by 1

        @param atoms: indexes of atoms that may have hydrogen attached
        @type atoms : list
        @return: indexes of attached hydrogen atoms
        @rtype : list
        """

        hydrogens = set()
        
        for j in atoms:
            atom = self.molecule.atoms[j].OBAtom
            for neighbor in pybel.ob.OBAtomAtomIter(atom):
                bond = atom.GetBond(neighbor)
                a_i = bond.GetBeginAtomIdx() - 1
                a_n = bond.GetBeginAtom().GetAtomicNum()
                b_i = bond.GetEndAtomIdx() - 1
                b_n = bond.GetEndAtom().GetAtomicNum()

                if a_n == 1:
                    hydrogens.add(a_i)
                if b_n == 1:
                    hydrogens.add(b_i)

        h = sorted(list(hydrogens))
        return h


class Geotool(object):
    def make_fragment(self, representation, kind="smiles"):
        """Take a linear representation of a molecule and add a title with
        IUPAC and SMILES designations. Convert a linear molecular specification
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

    def make_system(self, items, kind="smiles"):
        fragments = []
        if type(items) == str:
            items = [items]

        for item in items:
            fragment = self.make_fragment(item, kind=kind)
            fragments.append(fragment)

        s = System(fragments)
        return s

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
