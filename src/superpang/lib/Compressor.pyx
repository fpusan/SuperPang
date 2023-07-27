#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: cpow=True

from cpython.array cimport array
import numpy as np
cimport numpy as np
DTYPE = np.uint8
ctypedef np.uint8_t DTYPE_t

cdef class Compressor:
    """
    Hash a nucleotide sequence into a compressed form
    Each tetranucleotide will be stored as a single byte
    When the sequence length is not divisible by 4, each
     remaining base will be stored as a byte
    """

    cdef Py_ssize_t ksize
    cdef Py_ssize_t com
    cdef Py_ssize_t rem
    cdef Py_ssize_t clength

    def __init__(self, ksize):
        self.ksize = ksize
        self.com = self.ksize // 4
        self.rem = self.ksize % 4
        self.clength = self.com + self.rem
        
    def compress(self, const unsigned char [:] seq):
        cdef Py_ssize_t i, j
        cdef const unsigned char [:] tetra
        cdef DTYPE_t [1000] cseq
        cdef int ctetra
        for i in range(self.com):
            # Each tetranucleotide is stored as an 8-bit binary
            # A: 00, C: 01, G: 10, T: 11 
            tetra = seq[i*4:(i*4)+4]
            ctetra = 0
            for j in range(3,-1,-1):
                if tetra[3-j] == b'A':
                    pass
                elif tetra[3-j] == b'C':
                    ctetra += 10**(2*j)
                elif tetra[3-j] == b'G':
                    ctetra += 10**((2*j)+1)
                elif tetra[3-j] == b'T':
                    ctetra += 10**((2*j)+1)
                    ctetra += 10**(2*j)
            # Transform the 8-bit into a 1-byte decimal
            cseq[i] = self.change_base(ctetra, 2, 10)
        # For the remainder, just store each base as a byte
        for i in range(self.rem):
            cseq[self.com+i] = seq[4*self.com+i] # this directly casts the char into into DTYPE_t

        return cseq[0:self.clength]


    def decompress(self, const unsigned char[:] cseq):
        cdef Py_ssize_t i, j
        cdef int dec
        cdef int ctetra
        seq = array('u', '0'*self.ksize)
        for i in range(self.com):
            dec = cseq[i]
            # Get the 8-bit binary from the decimal byte
            ctetra = self.change_base(dec, 10, 2)
            # Reconstruc the tetranucleotide
            for j in range(3,-1,-1):
                if ctetra // 10**((2*j)+1):
                    if ( ctetra - 10**((2*j)+1) ) // 10**((2*j)):
                        seq[(i*4)+(3-j)] = 'T'
                        ctetra -= 10**((2*j)+1)
                        ctetra -= 10**(2*j)
                    else:
                        seq[(i*4)+(3-j)] = 'G'
                        ctetra -= 10**((2*j)+1)
                else:
                    if ctetra // 10**((2*j)):
                        seq[(i*4)+(3-j)] = 'C'
                        ctetra -= 10**(2*j)
                    else:
                        seq[(i*4)+(3-j)] = 'A'

        # For the remainder, just transform the byte into a character
        for i in range(self.rem):
            dec = cseq[self.com+i]
            seq[4*self.com+i] = chr(dec)
        
        return seq.tounicode()

    
    cdef inline int change_base(self, int num, int base0, int base1):
        # This expects 8-bit binary numbers and decimals no larger than 255!
        cdef Py_ssize_t i
        cdef int res = 0
        for i in range(7,-1,-1):
            if num // base1**i:
                num -= base1**i
                res += base0**i
        return res
