#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

from libc.stdlib cimport malloc, strtol
from cpython.array cimport array, clone

cdef Py_ssize_t MAX_RC_LEN = 50000000
cdef char *seq_dest = <char *>malloc(MAX_RC_LEN + 1)
seq_dest[MAX_RC_LEN] = b'\0'
cdef char *seq_destC = <char *>malloc(MAX_RC_LEN + 1)
seq_destC[MAX_RC_LEN] = b'\0'


cdef char *basemap = [ b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',   b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'T', b'\0', b'G', b'\0', b'\0', b'\0', b'C', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'N', b'\0', b'\0', b'\0', b'\0', b'\0', b'A', b'A', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',  b't', b'\0',  b'g', b'\0', b'\0', b'\0',  b'c', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',  b'a',  b'a' ]

# see Cython implementation (v2) from https://bioinformatics.stackexchange.com/questions/3583/what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho
cpdef reverse_complement(str seq):
    cdef bytes py_bytes = seq.encode('UTF-8')
    cdef char *seq_src = py_bytes
    cdef Py_ssize_t seq_len = len(py_bytes)
    if seq_len > MAX_RC_LEN:
        raise Exception('Too long input sequence, recompile!')
    cdef Py_ssize_t i = 0
    for i in range(seq_len):
        seq_dest[seq_len - i - 1] = basemap[<int>seq_src[i]]
    return seq_dest[:seq_len].decode('UTF-8')



cpdef parse_cigar(str cigar):
    cdef bytes py_bytes = cigar.encode('UTF-8')
    cdef Py_ssize_t cLen = len(py_bytes)
    cdef char *cigar_src = py_bytes
    cdef char c
    cdef char [10] buf = [0]*10 # Explicitly initialize to 0.
                                # Using conda compilers (but not base ubuntu compilers, even when versions were equal)
                                #  caused trailing characters being passed to strtol in some instances (buf seemed to have more than 10 elements)
                                #  so if those characters could be parsed as a number, the result from strtol was wront. Initializing to 0 from the starts apparently fixes this
    cdef Py_ssize_t i, j, bsize = 0, nOps = 0, idlen = 0
    cdef float mlen = 0
    cdef long [:] Larray
    cdef char [:] oparray
    cdef double iden
    Larray  = clone(array('l'), 100000, False)
    oparray = clone(array('b'), 100000, False)
    for i in range(cLen):
        c = cigar_src[i]
        if c >= 48 and c <= 57:
            buf[bsize] = c
            bsize += 1
        else:
            for j in range(bsize): # shift to the end
                buf[9-j] = buf[bsize-j-1]
            for j in range(10-bsize):
                buf[j] = 48 # "0"
            bsize = 0
            Larray[nOps] = strtol(buf, NULL, 10)
            oparray[nOps] = c
            if c == 61: # "="
                mlen  += Larray[nOps]
                idlen += Larray[nOps]
            elif c == 88: # "X":
                idlen += Larray[nOps]
            nOps += 1
    if idlen > 0:
        iden = mlen/idlen
    else:
        iden = 0
    return Larray[:nOps], oparray[:nOps], idlen, iden


cpdef correct_query(str query, Py_ssize_t queryStart, Py_ssize_t queryEnd, str target, Py_ssize_t targetStart, Py_ssize_t targetEnd, long [:] cigLengths, char[:] cigOps, bint isRC, int mismatch_size_threshold, int indel_size_threshold):
  
    query = query if not isRC else reverse_complement(query)

    cdef bytes q_bytes = query.encode('UTF-8')
    cdef bytes t_bytes = target.encode('UTF-8')
    cdef char *query_src  = q_bytes
    cdef char *target_src = t_bytes

    cdef Py_ssize_t nOps = cigLengths.size

    cdef Py_ssize_t queryPos, targetPos, resPos
    cdef Py_ssize_t i, j, L
    cdef char op

    queryPos = queryStart if not isRC else len(query) - queryEnd
    targetPos = targetStart
    resPos = 0

    for i in range(nOps):
        L  = cigLengths[i]
        op = cigOps[i]
        if op == 61:
            for j in range(L):
                if target_src[targetPos+j] != query_src[queryPos+j]:
                    raise ValueError
        if op == 61 or op == 88: # "=" or "X"
            if op == 61 or L <= mismatch_size_threshold:
                for j in range(L):
                    seq_destC[resPos+j] = target_src[targetPos+j]
                    assert target_src[targetPos+j] != 0
                resPos += L
            else:
                for j in range(L):
                    seq_destC[resPos+j] = query_src[queryPos+j]
                resPos += L
            queryPos += L
            targetPos += L
        elif op == 73: # "I"
            if L > indel_size_threshold:
                for j in range(L):
                    seq_destC[resPos+j] = query_src[queryPos+j]
                resPos += L
            queryPos += L

        elif op == 68: #"D"
            if L <= indel_size_threshold:
                for j in range(L):
                    seq_destC[resPos+j] = target_src[targetPos+j]
                resPos += L
            targetPos += L
        else:
            assert False

    #for j in range(resPos):
    #    assert seq_destC[j] != 0

    correction = seq_destC[:resPos].decode('UTF-8')
    
    if isRC:
        correction = reverse_complement(correction)

    return correction

