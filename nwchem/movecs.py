# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import sys
import struct
import argparse
import os
try:
    import cPickle as pickle
except ImportError:
    import pickle

#http://svn.assembla.com/svn/Zori/branches/DZori/python/nw_movecs.py
#nw_movecs.py:
#  extract mo coefficients from nwchem's .movecs file
#  a python port of Tim Heaton-Burgess' nwchem_movecs.F
#
#  nw_read_movecs returns a list of mo coefficient matrices: mo_set
#    C_alpha=mo_set[0];  C_beta=mo_set[1];  C[i][mu]
#
#  known bugs:
#   general: assume native endian-ness
#   frec(): portability: size of head and tail may vary with platform
#Brian Austin (May 5, 2008)

#comments above remain intact from original code
#many bug fixes/enhancements by Matt Ernst (April 26, 2014)

class MOVecsReader(object):
    """You can't use this class directly, only its subclasses, because it
    has no information about integer size or endianness.
    """
    
    def __init__(self, movecs_file_name, default_energy=0.0, default_enrep=0.0):
        """Initialize reader with the input file name and (optionally)
        default total energy and default nuclear repulsion energy.

        Before release 4.1 NWChem did not store these energies in movecs, so
        if writing a modern movecs file from old movecs input the user may
        need to supply the energy values manually, e.g. by copying them from
        the calculation log file. The default energies are used only
        if they cannot be found in the input files.

        @param movecs_file_name: name of molecular orbital vector file to read
        @type movecs_file_name : str
        @param default_energy: total energy in Hartrees
        @type default_energy : float
        @param default_enrep: nuclear repulsion energy in Hartrees
        @type default_enrep : float
        """
        
        self.fname = movecs_file_name
        self.default_energy = default_energy
        self.default_enrep = default_enrep

        #just test that the file can be opened, in case user gave bad
        #file name
        with open(self.fname, 'rb') as f:
            pass

    def _init_formats(self):
        """Build class-specific struct format strings for reading head,
        tail, and contents from Fortran records based on class
        settings of machine integer sizes and endianness.
        """

        try:
            if self._initialized:
                return
        except AttributeError:
            self._initialized = True
        
        if self.endianness == 'little':
            self._end = '<'

        elif self.endianness == 'big':
            self._end = '>'

        else:
            raise ValueError("Bad endianness {0}".format(self.endianness))

        #head and tail sizes expected as signed short,
        #signed long, or signed long long
        size_to_code = {2 : 'h',
                        4 : 'i',
                        8 : 'q'}

        self.head_fmt = self._end + size_to_code[self.head_size]
        self.tail_fmt = self._end + size_to_code[self.tail_size]
        self.int_fmt = size_to_code[self.int_size]

    def get_head(self, open_file):
        """Get the head of a Fortran record, basically a count of how
        many bytes are in the actual record data.

        @param open_file: open file to read a record head from
        @type open_file : file
        @return: number of bytes in record
        @rtype : int
        """
        
        head = open_file.read(self.head_size)
        head_unpacked = struct.unpack(self.head_fmt, head)
        return head_unpacked
        
    def get_tail(self, open_file):
        """Get the tail of a Fortran record. The actual value is unused
        but the read is necessary to advance the file. The tail appears
        to be a duplicate of the record head. Maybe to allow easy reading
        of the file in reverse order?

        @param open_file: open file to read a record head from
        @type open_file : file
        @return: number of bytes in record
        @rtype : int
        """
        
        tail = open_file.read(self.tail_size)
        tail_unpacked = struct.unpack(self.tail_fmt, tail)
        return tail_unpacked

    def get_record(self, open_file):
        """Get record data from a Fortran record within a .movecs file,
        discarding record tail.

        @param open_file: open file to read a record from
        @type open_file : file
        @return: record data as a raw byte string
        @rtype : str
        """
        head = self.get_head(open_file)

        #number of bytes to read from record came from the record head
        rec = open_file.read(head[0])

        #don't do anything with the record tail, just read it to advance
        #position in the file
        tail = self.get_tail(open_file)

        return rec

    def get_ints(self, open_file):
        """Read a record that consists of one or more integer from an
        open file. If there is only a single record, return as int
        instead of list.

        @param open_file: open file to read a record from
        @type open_file : file
        @return: an integer
        @rtype : int | list
        """
        
        raw_ints = self.get_record(open_file)
        num_ints = len(raw_ints) / self.int_size
        #the struct format string: endianness flag followed by one
        #signed int symbol for each int to be decoded from the raw data
        struct_fmt = self._end + (self.int_fmt * num_ints)
        unpacked = struct.unpack(struct_fmt, raw_ints)

        if len(unpacked) == 1:
            return unpacked[0]
        else:
            return unpacked

    def get_doubles(self, open_file):
        """Fetch data from a record that represents contiguous
        double precision floating point numbers. Convert the
        numbers to Python's float type (which is really a double) and
        return them in a list.

        @param open_file: open file to read a record from
        @type open_file : file
        @return: floating point numbers
        @rtype : list
        """
        double_size = 8
        
        raw_doubles = self.get_record(open_file)
        num_doubles = len(raw_doubles) / double_size

        #the struct format string: endianness flag followed by one 'd'
        #for each double to be decoded from the raw data
        struct_fmt = self._end + ('d' * num_doubles)
        unpacked = struct.unpack(struct_fmt, raw_doubles)
        return list(unpacked)

    def read(self):
        """Read the contents of .movecs file and perform any necessary
        platform translations.

        @return: structured, translated file contents
        @rtype : dict
        """
        
        self._init_formats()

        rdata = {'metadata' : {}}
        rdata['metadata']['reader_name'] = self.__class__.__name__
        rdata['metadata']['movecs_name'] = self.fname
        fsize = os.path.getsize(self.fname)
        
        with open(self.fname, 'rb') as f:
            ##c would be performed by movecs_read_header

            ##read(unitno) ! convergence info
            conv = self.get_record(f)

            #'conv' packs together a bunch of strings
            #basissum 32
            #geomsum 32
            #bqsum 32
            #scftype 20
            #date 26
            #the *sum values are md5 checksums, bqsum currently appears unused
            rdata['basissum'] = conv[:32]
            rdata['geomsum'] = conv[32:64]
            rdata['bqsum'] = conv[64:96]
            rdata['scftype'] = conv[96:116]
            rdata['date'] = conv[116:]

            ##read(unitno) ! scf type
            scf_type = self.get_record(f)
            rdata['scf_type'] = scf_type

            ##read(unitno) ! length of title
            ##read(unitno) ! title(1:lentit)
            title_len = self.get_ints(f)
            title = self.get_record(f)
            rdata['title'] = title

            ##read(unitno) ! length of basis name
            ##read(unitno) ! basis_name(1:lenbas)
            basis_name_len = self.get_ints(f)
            basis_name = self.get_record(f)
            rdata['basis_name'] = basis_name

            ##read(unitno) nsets ! = 1 for RHF/KS and 2 for UHF/KS
            Nmo_sets = self.get_ints(f)
            if Nmo_sets == 1:
                pNmo_sets = "Nmo_sets {0} : Restricted RHF/KS".format(Nmo_sets)
            elif Nmo_sets == 2:
                pNmo_sets = "Nmo_sets {0} : Unrestricted UHF/KS".format(Nmo_sets)
            else:
                pNmo_sets = "Nmo_sets {0} : unknown".format(Nmo_sets)

            rdata['wavefun_restrictions'] = pNmo_sets

            ##read(unitno) nbf   ! cardinality of orbital basis 
            ##write(*,*) 'nbf            ', nbf
            Nao = self.get_ints(f)
            rdata['Nbas_ao'] = Nao

            if Nmo_sets == 1:
                ##read(unitno) (nmo(i),i=1,nsets) ! nmo<nbf if linear dependencies
                Nmo = [self.get_ints(f)]
                rdata["Nbas_mo"] = Nmo[0]

            elif Nmo_sets == 2:
                ##read(unitno) (nmo(i),i=1,nsets) ! nmo<nbf if linear dependencies
                ##write(*,*) 'nmo alpha beta ', nmo(1),nmo(2)
                Nmo = list(self.get_ints(f))
                rdata["Nbas_mo: alpha"] = Nmo[0]
                rdata["Nbas_mo: beta"] = Nmo[1]

            mo_set = []
            ##do iset=1,nsets
            for iset in range(Nmo_sets):

                ##read(unitno) (occ(j,iset),j=1,nbf)   ! orbital occupancy
                orbN = self.get_doubles(f)

                ##read(unitno) (evals(j,iset),j=1,nbf) ! orbital energies
                orbE = self.get_doubles(f)

                ##! now MO coefficents
                ##do i=1,nmo(iset)
                ##  !orbitals form columns in movecs
                ##  read(unitno) ( movecs(j,i,iset),j=1,nmo(iset))
                ##enddo
                orbCv = []
                for j in range(Nao):
                    Cv = self.get_doubles(f)
                    orbCv += Cv
                orbC = []
                for k in range(Nmo[iset]):
                    orbC.append(orbCv[k * Nao : (k + 1) * Nao])

                mo_set.append([orbN, orbE, orbC])

            #If there are still bytes left in the file, process them.
            #Before 2001/12/28, NWChem release 4.1, nothing more stored at end.
            #After 2001/12/28, release 4.1, total energy stored.
            #After 2002/02/19, release 4.1, also nuclear repulsion energy.
            
            #If translating from old movecs, set missing values using
            #assigned defaults
            if f.tell() < fsize:
                energy, enrep = self.get_doubles(f)
                rdata['effective_nuclear_repulsion_energy'] = enrep
                rdata['total_energy'] = energy

            else:
                rdata['total_energy'] = self.default_energy
                rdata['effective_nuclear_repulsion_energy'] = self.default_enrep

            #enddo

            ##close(unitno)

        rdata['mo_set'] = mo_set

        return rdata

class MOVecsWriter(object):
    """You can't use this class directly, only its subclasses.

    Each writer is based on a MOVecsReader subclass for specific machine type.
    """

    def __init__(self, movecs_file_name):
        """Prepare writer for molecular orbital vectors.

        @param movecs_file_name: name of molecular orbital vector file to write
        @type movecs_file_name : str
        """
        
        self.fname = movecs_file_name
        
        #just test that the file can be opened for writing, in case user gave
        #bad file name or doesn't have permissions
        with open(self.fname, 'wb') as f:
            pass

        os.unlink(self.fname)

    def init_formats(self):
        """Build class-specific struct format strings for writing head,
        tail, and records into Fortran records based on class
        settings of machine integer sizes and endianness.
        """
        
        if self.reader.endianness == 'little':
            self._end = '<'

        elif self.reader.endianness == 'big':
            self._end = '>'

        #head and tail sizes expected as signed short,
        #signed long, or signed long long
        size_to_code = {2 : 'h',
                        4 : 'i',
                        8 : 'q'}

        self.head_fmt = self._end + size_to_code[self.reader.head_size]
        self.tail_fmt = self._end + size_to_code[self.reader.tail_size]
        self.int_fmt = size_to_code[self.reader.int_size]

    def put_head(self, open_file, count):
        """Put a head on a Fortran record, basically a count of how
        many bytes are in the actual record data.

        @param open_file: open file to write a record head to
        @type open_file : file
        @param count: number of bytes in record
        @type count : int
        """

        packed = struct.pack(self.head_fmt, count)
        open_file.write(packed)

    def put_record(self, open_file, packed):
        """Put a record into the open output file.

        @param open_file: open file to read a record head from
        @type open_file : file
        @param packed: a byte-string of record data packed by struct module
        @type packed : str
        """
        
        head_size = len(packed)

        self.put_head(open_file, head_size)
        open_file.write(packed)

        #tail is just head repeated after record data
        self.put_head(open_file, head_size)

    def put_string(self, open_file, s):
        """Write string into a record.

        @param open_file: open file to write a record to
        @type open_file : file
        @param s: string data
        @type s : str
        """

        fmt = 'c' * len(s)
        packed = struct.pack(fmt, *tuple(s))
        self.put_record(open_file, packed)

    def put_ints(self, open_file, ints):
        """Write one or more integers into a contiguous record.

        @param open_file: open file to write a record to
        @type open_file : file
        @param ints: string data
        @type ints : list
        """

        fmt = self._end + self.int_fmt * len(ints)
        packed = struct.pack(fmt, *tuple(ints))
        self.put_record(open_file, packed)

    def put_doubles(self, open_file, doubles):
        """Write one or more double precision floating point numbers into a
        contiguous record.

        @param open_file: open file to write a record to
        @type open_file : file
        @param doubles: double precision floating point values
        @type doubles : list
        """

        fmt = self._end + 'd' * len(doubles)
        packed = struct.pack(fmt, *tuple(doubles))
        self.put_record(open_file, packed)

    def write(self, data):
        """Write formatted data interpreted by a MOVecsReader into an NWChem
        .movecs file.

        @param data: formatted platform-independent data
        @type data : dict
        """

        with open(self.fname, 'wb') as f:
            #write convergence section data
            conv_elements = [data['basissum'], data['geomsum'], data['bqsum'],
                             data['scftype'], data['date']]
            conv_string = ''.join(conv_elements)
            self.put_string(f, conv_string)
            
            #write scf_type
            self.put_string(f, data['scf_type'])

            #write title with title length (which is actually fixed)
            self.put_ints(f, [len(data['title'])])
            self.put_string(f, data['title'])

            #write basis_name with title length (which is actually fixed)
            self.put_ints(f, [len(data['basis_name'])])
            self.put_string(f, data['basis_name'])

            #number of molecular orbital sets and their
            #respective orbital counts -- two sets for UHF, else one set
            if 'Unrestricted' in data['wavefun_restrictions']:
                Nmo_sets = 2
                Nmo = [data['Nbas_mo: alpha'],
                       data['Nbas_mo: beta']]

            else:
                Nmo_sets = 1
                Nmo = [data['Nbas_mo']]

            #write number of molecular orbital sets
            self.put_ints(f, [Nmo_sets])

            #write number of atomic orbitals
            Nao = data['Nbas_ao']
            self.put_ints(f, [Nao])

            #write number of molecular orbitals for each orbital set
            self.put_ints(f, Nmo)

            for iset in range(Nmo_sets):
                orbN, orbE, orbC = data['mo_set'][iset]

                #write orbital occupancy
                self.put_doubles(f, orbN)

                #write orbital energy
                self.put_doubles(f, orbE)

                #write MO coefficients
                for j in range(Nao):
                    Cv = orbC[j]
                    self.put_doubles(f, Cv)

            #write total energy and effective nuclear repulsion energy
            self.put_doubles(f, [data['total_energy'],
                                 data['effective_nuclear_repulsion_energy']])


class MOVecsReaderLittle64(MOVecsReader):
    """This is for reading molecular orbital vector files generated by NWChem
    on 64 bit little endian systems such as x86 processors running 64 bit
    Linux."""

    int_size = 8
    head_size = 4
    tail_size = 4
    endianness = 'little'
    
    def __init__(self, fname):
        MOVecsReader.__init__(self, fname)
        

class MOVecsReaderLittle32(MOVecsReader):
    """This is for files generated on 32 bit little endian systems
    such as x86 processors running 32 bit Linux."""

    int_size = 4
    head_size = 4
    tail_size = 4
    endianness = 'little'
    
    def __init__(self, fname):
        MOVecsReader.__init__(self, fname)
        

class MOVecsReaderBig32(MOVecsReader):
    """This is for files generated on 32 bit big endian systems
    such as MIPS processors running 32 bit Irix."""

    int_size = 4
    head_size = 4
    tail_size = 4
    endianness = 'big'
    
    def __init__(self, fname):
        MOVecsReader.__init__(self, fname)


class MOVecsReaderBig64(MOVecsReader):
    """This is for files generated on 64 bit big endian systems
    such as MIPS processors running 64 bit Irix."""

    int_size = 8
    head_size = 4
    tail_size = 4
    endianness = 'big'
    
    def __init__(self, fname):
        MOVecsReader.__init__(self, fname)
        

class MOVecsWriterLittle64(MOVecsWriter):
    """This is for writing molecular orbital vector files to match NWChem
    movecs format on 64 bit little endian systems such as x86 processors
    running 64 bit Linux."""

    def __init__(self, fname):
        MOVecsWriter.__init__(self, fname)
        self.reader = MOVecsReaderLittle64
        self.init_formats()

class MOVecsWriterLittle32(MOVecsWriter):
    """This is for writing molecular orbital vector files to match NWChem
    movecs format on 32 bit little endian systems such as x86 processors
    running 32 bit Linux."""

    def __init__(self, fname):
        MOVecsWriter.__init__(self, fname)
        self.reader = MOVecsReaderLittle32
        self.init_formats()

class MOVecsWriterBig64(MOVecsWriter):
    """This is for writing molecular orbital vector files to match NWChem
    expected input on 64 bit big endian systems such as MIPS processors
    running 64 bit Irix."""

    def __init__(self, fname):
        MOVecsWriter.__init__(self, fname)
        self.reader = MOVecsReaderBig64
        self.init_formats()

class MOVecsWriterBig32(MOVecsWriter):
    """This is for writing molecular orbital vector files to match NWChem
    movecs format on 32 bit big endian systems such as MIPS processors
    running 32 bit Irix."""

    def __init__(self, fname):
        MOVecsWriter.__init__(self, fname)
        self.reader = MOVecsReaderBig32
        self.init_formats()

class PickleReader(object):
    """This is for reading previously translated, pickled files. It's provided
    to avoid special casing readers and writers when translating files.
    """

    def __init__(self, fname):
        self.fname = fname

    def read(self):
        with open(self.fname, 'rb') as f:
            data = pickle.load(f)

        return data

class PickleWriter(object):
    """This is for serializing previously translated vectors as a platform
    independent pickle file. It's provided to avoid special casing readers
    and writers when translating files.
    """

    def __init__(self, fname):
        self.fname = fname

    def write(self, data):
        with open(self.fname, 'wb') as f:
            pickle.dump(data, f)

def native_platform():
    """Guess which kind of platform this script is currently running on.
    It could guess wrong if you're doing something tricky like running a
    32 bit Python interpreter under a 64 bit OS.

    @return: a name representing platform endianness and bitness, like 'Big64'
    @rtype : str
    """

    endianness = sys.byteorder.title()
    if sys.maxsize == 2 ** 63 - 1:
        bitness = '64'

    elif sys.maxsize == 2 ** 31 - 1:
        bitness = '32'

    else:
        raise ValueError("Native integer doesn't appear to be 32 or 64 bits -- this shouldn't happen")

    pname = endianness + bitness
    return pname

class GuessReader(object):
    def __init__(self, fname):
        self.fname = fname

    def read(self):    
        """Throw all readers at a file, hope that exactly one sticks. Send
        back translated data if there is no ambiguity, else raise an error.

        @return: translated data
        @rtype : dict
        """

        readers = [MOVecsReaderLittle64, MOVecsReaderLittle32,
                   MOVecsReaderBig32, MOVecsReaderBig64, PickleReader]
        readers = [readers[3]]

        results = []
        for reader in readers:
            rdata = None
            R = reader(self.fname)
            try:
                rdata = R.read()
            except Exception, e:
                pass

            results.append(rdata)

        ok_readers = len(results) - results.count(None)
        names = [x['metadata']['reader_name'] for x in results if x is not None]

        if ok_readers > 1:
            err = "Ambiguous format -- multiple readers worked: {0}".format(names)
            raise ValueError(err)

        elif ok_readers == 1:
            data = [x for x in results if x is not None][0]
            return data

        else:
            raise ValueError("Failed to read file")

if __name__ == '__main__':
    cmap = {'Little64' : (MOVecsReaderLittle64, MOVecsWriterLittle64),
            'Little32' : (MOVecsReaderLittle32, MOVecsWriterLittle32),
            'Big64' : (MOVecsReaderBig64, MOVecsWriterBig64),
            'Big32' : (MOVecsReaderBig32, MOVecsWriterBig32),
            'Pickle' : (PickleReader, PickleWriter)}
    
    ckeys = sorted(cmap.keys())

    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile', help="File to translate, either a .movecs file or a previously translated .pkl file.")
    parser.add_argument('outputfile', help="Name of translated file to write, either a .movecs file for NWChem or a pickled Python data structure file ending in .pkl")
    parser.add_argument('-r', '--reader', default='guess', help="Kind of reader to use on an input .movecs or pickle file, one of {0}. The default is to guess based on the file contents. MOVecs files created on an x86, Alpha, or Itanium system running Linux are little-endian. Those created on PowerPC, SPARC, or MIPS are big-endian. 32 or 64 bit integer size depends on CPU, OS, and compiler. Pickle files are platform-independent.".format(ckeys))
    parser.add_argument('-w', '--writer', default='guess', help="Kind of writer to use storing an output file, one of {0}. The default is to guess based on the output file extension and (possibly) the current execution platform.".format(ckeys))
    parser.add_argument('--enrep', default=0.0, help="Nuclear repulsion energy to assign given in atomic units (Hartree), if input .movecs file does not contain the nuclear repulsion energy. NWChem did not store this value in movecs files before release 4.1.")
    parser.add_argument('--energy', default=0.0, help="Default total energy to assign given in atomic units (Hartree), if input .movecs file does not contain the total energy. NWChem did not store this value in movecs files before release 4.1.")
    args = parser.parse_args()

    if args.inputfile == args.outputfile:
        sys.stderr.write("Can't write to the same file you are reading from\n")
        sys.exit(1)

    if args.reader == 'guess':
        reader = GuessReader

    else:
        try:
            reader = cmap[args.reader][0]
        except KeyError:
            sys.stderr.write("Invalid reader name\n")
            parser.print_help()
            sys.exit(1)

    if args.writer == 'guess':
        if args.outputfile.endswith('.pkl'):
            writer = PickleWriter

        else:
            platform = native_platform()
            writer = cmap[platform][1]
            
    else:
        try:
            writer = cmap[args.writer][1]
        except KeyError:
            sys.stderr.write("Invalid writer name\n")
            parser.print_help()
            sys.exit(1)

    try:
        data = reader(args.inputfile).read()
        W = writer(args.outputfile)
        W.write(data)
    except Exception, e:
        sys.stderr.write(str(e) + '\n')
