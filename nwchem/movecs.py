# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import sys
import pprint
import struct
import array
import os
try:
    import cPickle as pickle
except ImportError:
    import pickle

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

#header above remains intact from original code
#bug fixes/enhancements by Matt Ernst (April 26, 2014)

class MOVecsReader(object):
    """You can't use this class directly, only its subclasses, because it
    has no information about record integer size or endianness.
    """
    
    def __init__(self, movecs_file_name):
        self.fname = movecs_file_name

        #just test that the file can be opened, in case user gave bad
        #file name
        with open(self.fname, 'r') as f:
            pass

    def _init_formats(self):
        """Build class-specific struct format strings for reading head,
        tail, and integers from raw Fortran records based on class
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

        #head and tail sizes expected as unsigned short,
        #unsigned long, or unsigned long long
        size_to_code = {2 : 'H',
                        4 : 'I',
                        8 : 'Q'}

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
        #unsigned int symbol for each int to be decoded from the raw data
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
        self._init_formats()
        
        rdata = {'reader_name' : self.__class__.__name__}
        rdata['movecs_name'] = self.fname
        fsize = os.path.getsize(self.fname)
        
        with open(self.fname, 'rb') as f:
            ##c would be performed by movecs_read_header

            ##read(unitno) ! convergence info
            conv = self.get_record(f)
            rdata['conv'] = conv.strip()

            ##read(unitno) ! scf type
            scf_type = self.get_record(f)
            rdata['scf_type'] = scf_type.strip()

            ##read(unitno) ! length of title
            ##read(unitno) ! title(1:lentit)
            title_len = self.get_ints(f)
            title = self.get_record(f)
            rdata['title'] = title.strip()

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

            if( Nmo_sets == 1 ):
                ##read(unitno) (nmo(i),i=1,nsets) ! nmo<nbf if linear dependencies
                Nmo = [self.get_ints(f)]
                rdata["Nbas_mo"] = Nmo[0]

            elif( Nmo_sets == 2 ):
                ##read(unitno) (nmo(i),i=1,nsets) ! nmo<nbf if linear dependencies
                ##write(*,*) 'nmo alpha beta ', nmo(1),nmo(2)
                Nmo = list(self.get_ints(f))
                rdata["Nbas_mo: alpha"] = Nmo[0]
                rdata["Nbas_mo: beta"] = Nmo[1]

            mo_set = []
            ##do iset=1,nsets
            for iset in range( Nmo_sets ):

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
                for i in range(Nao):
                    Cv = self.get_doubles(f)
                    orbCv += Cv
                orbC = []
                for i in range(Nmo[iset]):
                    #orbC.append( orbCv[i::Nao] )
                    orbC.append(orbCv[i * Nao : (i + 1) * Nao])
                orbCv = []

                mo_set.append([orbN, orbE, orbC])


            #If there are still bytes left in the file, it contains
            #two energy values. Very old files do not contain these
            #values.
            if f.tell() < fsize:
                energy, enrep = self.get_doubles(f)
                rdata['effective_nuclear_repulsion_energy'] = enrep
                rdata['total_energy'] = energy

            #enddo

            ##close(unitno)

        rdata['mo_set'] = mo_set

        return rdata

class MOVecsReaderLittle64(MOVecsReader):
    """This is for molecular orbital vector files generated by NWChem
    on 64 bit little endian systems such as x86 processors running 64 bit
    Linux."""
    def __init__(self, fname):
        MOVecsReader.__init__(self, fname)
        self.int_size = 8
        self.head_size = 4
        self.tail_size = 4
        self.endianness = 'little'

class MOVecsReaderLittle32(MOVecsReader):
    """This is for files generated on 32 bit little endian systems
    such as x86 processors running 32 bit Linux."""
    def __init__(self, fname):
        MOVecsReader.__init__(self, fname)
        self.int_size = 4
        self.head_size = 4
        self.tail_size = 4
        self.endianness = 'little'

class MOVecsReaderBig32(MOVecsReader):
    """This is for files generated on 32 bit big endian systems
    such as MIPS processors running 32 bit Irix."""
    def __init__(self, fname):
        MOVecsReader.__init__(self, fname)
        self.int_size = 4
        self.head_size = 4
        self.tail_size = 4
        self.endianness = 'big'

class MOVecsReaderBig64(MOVecsReader):
    """This is for files generated on 64 bit big endian systems
    such as MIPS processors running 64 bit Irix."""
    def __init__(self, fname):
        MOVecsReader.__init__(self, fname)
        self.int_size = 8
        self.head_size = 4
        self.tail_size = 4
        self.endianness = 'big'


if __name__ == '__main__':
    #throw all readers at a file, hope that exactly one sticks
    readers = [MOVecsReaderLittle64, MOVecsReaderLittle32,
               MOVecsReaderBig32, MOVecsReaderBig64]

    results = []
    for reader in readers:
        data = None
        R = reader(sys.argv[1])
        try:
            data = R.read()
        except Exception, e:
            pass

        results.append(data)

    ok_readers = len(results) - results.count(None)
    names = [x['reader_name'] for x in results if x is not None]
        
    if ok_readers > 1:
        sys.stderr.write("Ambiguous format -- multiple readers worked: {0}\n".format(names))

    elif ok_readers == 1:
        data = [x for x in results if x is not None][0]
        #pprint.pprint(data)
        #pprint.pprint(data.keys())

        pickled_name = sys.argv[1] + '.pkl'
        with open(pickled_name, 'w') as outfile:
            pickle.dump(data, outfile)

    else:
        sys.stderr.write("Failed to read file\n")
        
    
