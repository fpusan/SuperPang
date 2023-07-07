#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

cimport cython
from cpython.array cimport array, clone
cdef extern from "Python.h":
     cdef Py_ssize_t PY_SSIZE_T_MAX
import numpy as np
cimport numpy as np
np.import_array()

DTYPE = np.uint32
ctypedef np.uint32_t DTYPE_t


cpdef compress_vertices(DTYPE_t [:] vs):
    """
    Transform a numpy array containing a list of vertices into a compressed form
      storing intervals of vertices that were present in the list, in the form
       [[int1_start, int1_end],
        [int2_start, int2_end],
        ...
        [intN_start, intN_end]]
    """
    cdef DTYPE_t [:,:] res
    cdef Py_ssize_t maxz = vs.shape[0]//2 + 3
    res = np.empty(shape = ( maxz, 2), dtype = DTYPE)
    cdef bint firstvertex = True
    cdef Py_ssize_t i,v
    cdef DTYPE_t last = 0
    cdef DTYPE_t cs = 0
    cdef Py_ssize_t z = 0
    np.asarray(vs).sort()
    for i in range(vs.shape[0]):
        v = vs[i]
        if v - last > 1:
            if not firstvertex:
                res[z][0] = cs
                res[z][1] = last
                z += 1
            cs = v
        firstvertex = False
        last = v
    res[z][0] = cs
    res[z][1] = last
    assert z < maxz
    return np.asarray(res[:z+1,])




cpdef vertex_overlap(DTYPE_t[:,:] x, DTYPE_t[:,:] y):
    """
    Return the number of overlapping vertices between two compressed vertex lists (seqPaths)
    """
    cdef Py_ssize_t rx = 0
    cdef Py_ssize_t ry = 0
    cdef Py_ssize_t ol = 0
    cdef DTYPE_t xs, xe, ys, ye
    while (rx < x.shape[0] and ry < y.shape[0]):
        xs = x[rx,0]
        xe = x[rx,1]
        ys = y[ry,0]
        ye = y[ry,1]
        if ys > xe:   # x is behind, advance rx
            rx += 1
        elif xs > ye: # y is behind, advance ry
            ry += 1
        else:         # there is some overlap
            if xs >= ys:     # x starts before
                if xe < ye:  #  ... and ends before too
                    ol += (xe - xs + 1)
                    rx += 1
                else:        #  ... and y ends before
                    ol += (ye - xs + 1)
                    ry += 1
            else:            # y starts before
                if ye < xe:  #  ... and ends before too
                    ol += (ye - ys + 1)
                    ry += 1
                else:        #  ... and x ends before
                    ol += (xe - ys + 1)
                    rx += 1
    return ol


cpdef isInSeqPath(DTYPE_t v, DTYPE_t[:,:] vs):
    """
    Return whether a vertex is present in a compressed vertex list (seqPath)
    """
    cdef Py_ssize_t rvs = 0
    cdef DTYPE_t vss, vse
    cdef bint is_present = False
    while (rvs < vs.shape[0]):
        vss = vs[rvs,0]
        vse = vs[rvs,1]
        if v >= vss and v <= vse:
            is_present = True
            break
        rvs += 1
    return is_present

