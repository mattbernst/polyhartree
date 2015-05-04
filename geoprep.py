# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import os
import string
import copy
from cinfony import pybel, webel
from sharedutilities import ELEMENTS

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
        self.explicit_atom_properties = {}

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

    def atom_properties(self, name):
        """Get named atom properties from explicit system properties or from
        the underlying fragments in the system.

        @return: per-atom properties across all atoms
        @rtype : list
        """

        fragment_props = []
        for f in self.fragments:
            try:
                p = f.atom_properties[name]
            except KeyError:
                p = [None] * len(f.atoms)
                
            fragment_props += p

        e_props = self.explicit_atom_properties.get(name,
                                                    [None] * len(self.atoms))
        properties = []

        #try to get properties from explicit system-level values first, then
        #attempt to fill in from fragment properties if explicit data absent
        for k in range(len(e_props)):
            e_p = e_props[k]
            if e_p is not None:
                properties.append(e_p)
            else:
                f_p = fragment_props[k]
                properties.append(f_p)

        return properties

    def set_properties(self, name, selection, properties):
        """Assign properties grouped by name to selected atoms. Any unselected
        atoms will get a None property.

        @param name: name of property group, e.g. "basis_name", "fukui_mu_n(+)"
        @type name : str
        @param selection: atom indices
        @type selection : list
        @param properties: any sequence of values, length same as selection
        @type properties : list
        """

        slen = len(selection)
        plen = len(properties)
        if slen != plen:
            raise ValueError("Got selection of {0} atoms but {1} properties: {2} {3}".format(slen, plen, selection, properties))

        #initialize property group with a list of None
        if name not in self.explicit_atom_properties:
            self.explicit_atom_properties[name] = [None] * len(self.atoms)

        for j in range(slen):
            k = selection[j]
            value = properties[j]
            self.explicit_atom_properties[name][k] = value

    def select(self, smarts, hydrogen="include"):
        """Select atoms matching a SMARTS pattern with different treatments
        for hydrogen atoms, and do it across all fragments in the system.

        @param smarts: a SMARTS pattern
        @type smarts : str
        @param hydrogen: "include", "exclude", or "only"
        @type hydrogen : str
        @return: indexes of atoms matching the selection
        @rtype : list
        """

        selected = []
        natoms_last = 0
        for f in self.fragments:
            s = f.select(smarts, hydrogen=hydrogen)
            s = [k + natoms_last for k in s]
            selected += s
            natoms_last += len(f.atoms)

        return selected

    @property
    def elements(self):
        """Get all elements incorporated in system.

        @param system: System
        @type system : geoprep.System
        @return: unique list of system elements, in order of atomic number
        @rtype : list
        """
        
        included_elements = set(self.atom_properties("symbols"))

        indexed = [(ELEMENTS.index(e), e) for e in included_elements]
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

        self.bridged_attrs = ['addh', 'atoms', 'calcdesc', 'calcfp', 'charge',
                              'conformers', 'data', 'dim', 'draw', 'energy',
                              'exactmass', 'formula', 'localopt', 'make3D',
                              'molwt', 'removeh', 'sssr', 'write']
        self.molecule = molecule
        self.atom_properties = {}
        
        #keep a "universal SMILES" representation of the fragment handy
        m = pybel.Molecule(molecule)
        self.smiles = m.write("smi", opt={"U" : True}).strip()

        self._bind_bridged_attrs()
        self.add_title()
        self.assign_elements()

    def __deepcopy__(self, memo):
        """Create a deep copy of the fragment. The main issue is that the
        underlying molecule is not pure Python and needs special handling
        to duplicate.

        @param memo: memoization dictionary
        @type memo : dict
        @return: fragment duplicate
        @rtype : Fragment
        """

        #create a new fragment and then replace the details using self
        result = Geotool().make_fragment("O")
        memo[id(self)] = result
        
        for k, v in self.__dict__.items():
            if k == "molecule":
                m = pybel.Molecule(self.molecule)
                #atom coordinates lose precision above, and need manual copy
                for j, entry in enumerate(self.geometry_list):
                    coords = entry[1:]
                    m.atoms[j].OBAtom.SetVector(*coords)
                result.molecule = m

            elif k in self.bridged_attrs:
                pass
                
            else:
                setattr(result, k, copy.deepcopy(v, memo))

        result._bind_bridged_attrs()
            
        return result

    def __eq__(self, other):
        """Equality comparison. To be equal, two fragments must compare equal
        across all key attributes.

        @param other: other fragment to compare to
        @type other : Fragment
        @return: True if fragments are equal, else False
        @rtype : bool
        """

        key_attrs = ["spin", "charge", "geometry_list", "atom_properties",
                     "title", "smiles"]
        
        for k in key_attrs:
            s = getattr(self, k, None)
            o = getattr(other, k, None)
            if s != o:
                return False

        return True

    def __ne__(self, other):
        """Inequality comparison. This is just the logical inverse of equality.
        """

        return not (self == other)

    def _bind_bridged_attrs(self):
        """Try to pass through attributes/methods of the underlying molecule
        as attributes/methods of the fragment.
        """

        for attr in self.bridged_attrs:
            try:
                p = getattr(self.molecule, attr)
                setattr(self, attr, p)
            except AttributeError:
                pass

    def add_title(self):
        """Generate a title for the fragment with an IUPAC name (if we can
        look one up) and a SMILES representation.

        TODO: Cache webel calls. Re-enable IUPAC name lookup.
        """
        
        #iupac = webel.Molecule(self.molecule).write('iupac').strip()
        iupac = "NOT A REAL IUPAC NAME"
        self.title = "{0} SMILES: {1}".format(iupac, self.smiles)

    def assign_elements(self):
        """Assign element symbols to all atoms as atom_properties["symbols"]
        """

        symbols = [ELEMENTS[atom.atomicnum - 1] for atom in self.atoms]
        self.set_properties("symbols", range(len(symbols)), symbols)

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

        This may make Mopac7 a little less finicky with linear molecules.
        """

        base = self.molecule.atoms[0].coords
        neg = [a * -1 for a in base]
        self.translate(neg)

    @property
    def geometry_list(self):
        """Provide fragment geometry as a simple list of lists, e.g.
        [["H", 0.0, 0.0, 0.0], ["H", 0.9, 0.0, 0.0]]

        @return: geometry list
        @rtype : list
        """

        symbols = self.atom_properties["symbols"]
        g = []
        for k in range(len(self.molecule.atoms)):
            coords = list(self.molecule.atoms[k].coords)
            entry = [symbols[k]] + coords
            g.append(entry)

        return g

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

        geolist = self.geometry_list
        for k in range(len(geolist)):
            coords = geolist[k][1:]
            translated = []
            for j in range(3):
                t = coords[j] + vec[j]
                translated.append(t)

            self.molecule.atoms[k].OBAtom.SetVector(*translated)

    def write_fragment(self, name=None, fmt=None, handle=None):
        """Write fragment contents to named file. If fmt is not provided
        explicitly it will be derived from the file name. If handle is not
        provided explicitly, it will be generated by opening an output file
        with the given name.

        @param name: optional name of file to write
        @type name : str
        @param fmt: if provided, use this format for output
        @type fmt : str
        @param handle: if provided, use this output file handle to write data
        @type handle : file
        """

        if fmt is None:
            fmt = name.rsplit(".", 1)[-1]

        data = self.write(fmt)
        if handle is None:
            with open(name, "wb") as handle:
                handle.write(data)
        else:
            handle.write(data)

    def select(self, smarts, hydrogen="include", flatten=True):
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
        @param flatten: if True, merge all selections
        @type flatten : bool
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

        if flatten:
            results = sum(results, [])
            
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

    def set_properties(self, name, selection, properties):
        """Assign properties grouped by name to selected atoms. Any unselected
        atoms will get a None property.

        @param name: name of property group, e.g. "basis_name", "fukui_mu_n(+)"
        @type name : str
        @param selection: atom indices
        @type selection : list
        @param properties: any sequence of values, length same as selection
        @type properties : list
        """

        slen = len(selection)
        plen = len(properties)
        if slen != plen:
            raise ValueError("Got selection of {0} atoms but {1} properties: {2} {3}".format(slen, plen, selection, properties))

        #initialize property group with a list of None
        if name not in self.atom_properties:
            self.atom_properties[name] = [None] * len(self.atoms)

        for j in range(slen):
            k = selection[j]
            value = properties[j]
            self.atom_properties[name][k] = value

    def set_basis_name_general(self, selection, mapfn, **kw):
        """Assign basis set names to selected atoms.

        @param selection: initial selection (all atoms if empty)
        @type selection : list
        @param mapfn: a method that produces names for selection
        @type mapfn : method
        @param **kw: optional keyword arguments for mapfn
        @type **kw : dict
        """
        
        name = "basis_name"
        values = []

        #no selection is equivalent to selecting everything
        if not selection:
            selection = range(len(self.atoms))

        atoms = self.atoms
        for index in selection:
            atom = atoms[index]
            v = mapfn(atom, **kw)
            values.append(v)

        self.set_properties(name, selection, values)

    def set_basis_name(self, basis_name, selection=[]):
        """Simple name based assignment: same basis name on every atom in
        selection.

        @param basis_name: basis set name for selection
        @type basis_name : str
        @param selection: initial selection (all atoms if empty)
        @type selection : list
        """

        def namer(atom, **kw):
            basis_name = kw["basis_name"]
            return basis_name

        kw = {"basis_name" : basis_name}
        self.set_basis_name_general(selection, namer, **kw)

class Geotool(object):
    def make_fragment(self, representation, fmt="smiles"):
        """Take a linear representation of a molecule and add a title with
        IUPAC and SMILES designations. Convert a linear molecular specification
        to a 3D form.

        @param representation: linear molecule encoding
        @type representation : str
        @param fmt: smiles, inchi, etc. (default smiles)
        @type fmt : str
        @return: 3D-form molecular fragment
        @rtype: Fragment
        """

        molecule = pybel.readstring(fmt, representation)
        
        molecule.make3D()
        fragment = Fragment(molecule)
        fragment.set_zero_to_origin()
        
        return fragment

    def make_system(self, items, fmt="smiles"):
        """Make a system out of one or more linear representations of
        molecules that will become fragments.

        @param items: one or more linear molecule representations
        @type items : str | list
        @return: a system containing one or more fragments
        @rtype : System
        """
        
        fragments = []
        if type(items) == str:
            items = [items]

        for item in items:
            fragment = self.make_fragment(item, fmt=fmt)
            fragments.append(fragment)

        s = System(fragments)
        return s

    def read_fragment(self, name=None, fmt=None, handle=None):
        """Read a molecular structure from a file. Guess at the file type from
        extension if caller does not supply explicit fmt. If file handle is
        provided, read data from it. Otherwise open file name for reading.

        @param name: file to open
        @type name : str
        @param fmt: optional OpenBabel format code e.g. "xyz"
        @type fmt : str
        @return: molecular fragment
        @rtype : Fragment
        """

        if not fmt:
            try:
                fmt = name.rsplit(".", 1)[-1]
            except (IndexError, AttributeError):
                msg = "No fmt given for {0} and unable to guess from file extension".format(repr(name))

        if handle is None:
            molecule = pybel.readfile(fmt, name).next()
        else:
            data = handle.read()
            molecule = pybel.readstring(fmt, data)
            
        fragment = Fragment(molecule)
        fragment.set_zero_to_origin()
        
        return fragment

    def mcs(self, fragments):
        """Find the maximum common substructure from a list of fragments.

        N.B.: Currently does not expose the many options provided by rdkit:
        http://www.rdkit.org/Python_Docs/rdkit.Chem.MCS-module.html

        Also, SMARTS match naturally includes heavy atoms only.

        @param fragments: two or more fragments containing common substructure
        @type fragments : list
        @return: maximum common substructure result
        @rtype : MCSResult
        """

        try:
            global MCS
            global rdk
            MCS
        except NameError:
            from rdkit.Chem import MCS
            from cinfony import rdk

        rf = [rdk.Molecule(f.molecule).Mol for f in fragments]
        cs = MCS.FindMCS(rf)
        
        return cs

    def align(self, reference, target):
        """Optimize target fragment alignment to reference. Return a copied
        fragment with aligned geometry and information about the fit.

        TODO (maybe): use maximum common substructure to also align different
        molecules, e.g. toluene against m-xylene

        @param reference: the fragment to align to
        @type reference : Fragment
        @param target: the fragment to be aligned
        @type target: Fragment
        @return: aligned fragment and fit data
        @rtype : dict
        """

        copied = copy.deepcopy(target)
        aligner = pybel.ob.OBAlign(False, False)
        aligner.SetRefMol(reference.molecule.OBMol) 
        aligner.SetTargetMol(copied.molecule.OBMol) 
        aligner.Align() 
        aligner.UpdateCoords(copied.molecule.OBMol)

        r = {"rmsd" : aligner.GetRMSD(), "fragment" : copied}
        return r
