import numpy as np
cimport numpy as np
import time
from time import clock
from spinwaves.spinwavecalc.readfiles import atom, readFiles
from periodictable import elements

cimport cMcQueenyAlg as calg


"""
HOW THIS WORKS

This code was adapted from the code in mcQueeny_alg.py.  As of right now (3/24/11)
It does aexactly the same thing.  I rewrote the majority of the computations in
C and used cython to call the C functions.  The C functions that I wrote are
contained in the file, McQueeny_Alg.c and the "public" functions (those that I
call through cython), are declared in McQueeny_Alg.h.  Cython then requires that
I have something like a cython header file (.pxd), which makes these C functions
available in my cython code.  I named this file cMcQueenyAlg.pxs and it's
contents is almost identical to McQueeny_Alg.h, except for the Atom struct.  I
could have left that identical as well, but since the details are not important
to my implementation on the cython side, I can just put "pass" for the internal
details of the struct.  I can then import cMcQueenyAlg like I would a python
file, and use the methods and type defined in it.

This cython file is then comipled using distutils(There are other ways too) by
entering "python setupy.py build_ext --inplace".  This creates a C file with the
same name as the cython .pyx file, "mcQueeny_C_alg.c" and compiles it to a .pyd
file, "mcQueeny_C_alg.pyd".  This can then be imported like any other python
file and the methods in it can be called.  mcQueeny_case(...) is called in from
this file (well the compiled .pyd version of this file) the exact same way it is
called from "McQueeny_Alg.py", just by changing the name of the imported file.

I left the outer loop which loops through Q points in cython.  I also left the
eigenvalue and eigenvector calculation in cython, so that I can use numpy
rather than finding a C library for eignevectors or C binding to LAPACK.
I wrote various methods at the top of this file for converting between Python
types and the C types I defined.

Tom Sarvey 3/24/11
"""

cdef int* intListToCArray(list):
    cdef int* array = calg.allocateIntArray(len(list))
    for i in range(len(list)):
        array[i] = list[i]
    return array

cdef calg.Real* dblListToCArray(list):
    cdef calg.Real* array = calg.allocateRealArray(len(list))
    for i in range(len(list)):
        array[i] = list[i]
    return array

#takes a list of python SimpleAtoms and returns a C array of calg.Atom
cdef calg.Atom* createCAtomList(list, jnums, jmats, spinAxis):
    cdef calg.Atom *clist = calg.createAtomList(len(list))
    #cdef calg.Real *dblJs
    for i in range(len(list)):
        dblJs = []
        for j in range(len(list[i].interactions)):
            dblJs.append(get_Jij(jnums, jmats, list[i].interactions[j])[1][1])
        calg.addToList(clist, i, len(list[i].neighbors), intListToCArray(list[i].neighbors), dblListToCArray(dblJs), list[i].spin[spinAxis], list[i].pos[0], list[i].pos[1], list[i].pos[2])
        #dblJs = calg.allocateRealArray(len(list[i].interactions))
        #for j in range(len(list[i].interactions)):
        #    print "here"
        #    dblJs[j] = get_Jij(jnums, jmats, list[i].interactions[j])
        #    print "here1"
        #calg.addToList(clist, i, len(list[i].neighbors), intListToCArray(list[i].neighbors), dblJs, list[i].spin[spinAxis], list[i].pos[0], list[i].pos[1], list[i].pos[2])
    return clist


#convert a calg.SquareMatrix to a 2D numpy array (complex)
cdef np.ndarray[np.complex128_t, cast=True] sqr_mat_to_np(calg.SquareMatrix *sqMat):
    cdef np.ndarray[np.complex128_t, ndim=2] a = np.zeros(((sqMat[0]).dim, (sqMat[0]).dim), dtype=np.complex128)
    cdef int i,j
    for i in range((sqMat[0]).dim):
        for j in range((sqMat[0]).dim):
            #I need to typecast the unknown "Real" to double to appease Cython
            a[i][j] = <double>((sqMat[0]).mat[i][j].real) + 1j*(<double>((sqMat[0]).mat[i][j].imaginary))
    return a


#taken directly from mcQueeny_alg.py
def get_Jij(jnums, jmats, num):
    """Returns the 3x3 numpy array corresponding to the matrix number, num from an atoms list
    of interactions."""
    index = jnums.index(num)
    return np.float64(np.array(jmats[index])) #This would be a sympy Matrix if I didn't use the np wrapper

cdef void convert_evec(calg.ComplexNum **cEvec, evec):
    """Converts the numpy array of eigenvectors to a 2D C array of calg.ComplexNums (ComplexNum **)."""
    #cdef calg.ComplexNum **cEvec = calg.allocate_evec(len(evec))
    cdef int i,j
    for i in range(len(evec)):
        for j in range(len(evec)):
            cEvec[i][j].real = evec[i][j].real
            cEvec[i][j].imaginary = evec[i][j].imag


#Copy a python list,w (list of eigenvalues) to the calg.Real array cW.
#these should be real, so that is checked and an error is raised if there is
#a complex value.  They all have tiny complex components, so a somewhat
#arbitrary method is used to determine if the imaginary component is significant.
cdef void convert_w(calg.Real *cW, w):
    cdef int i
    cdef calg.Real avg
    #This isn't fast, but I think its a safer way of determining how small a
    #number has to be to be negligible than hardcoding a small constant
    for val in w:
        avg += np.abs(val)
    avg = avg/len(w)

    for val in w:
        if val.imag > (val.real*1e-8) and val.imag > avg/1000:
            print val
            raise Exception("Complex Eigenvalue Error!")
    #return dblListToCArray(w)
    for i in range(len(w)):
        cW[i] = w[i]

#only necessary for testing
##cdef np.ndarray[np.complex128_t, cast=True] c_evec_to_numpy(calg.ComplexNum **cEvec, int size):
##    cdef int i,j
##    cdef np.ndarray[np.complex128_t, ndim=2] evec = np.zeros((size, size), dtype=np.complex128)
##    for i in range(size):
##        for j in range(size):
##            evec[i][j] = cEvec[i][j].real + 1j*cEvec[i][j].imaginary
##    return evec


def mcQueeny_case(interactionfile, spinfile, tau_list, temperature,
                      direction={'kx':1,'ky':0,'kz':0}, hkl_interval=[1e-3,2*np.pi,100],
                      omega_interval=[0,5,100], eig_eps = 0.001):

    #benchmarking
    initial_time = time.clock()

    #read the data files
    atom_list, jnums, jmats, N_atoms_uc, numCutOff, magCellSize = readFiles(interactionfile,spinfile, rtn_cutOffInfo = True)

    #Convert the Python list of SimpleAtoms to a C array of calg.Atom
    cdef calg.Atom *clist = createCAtomList(atom_list, jnums, jmats, 2)

    #generate a few lists in Python (could be done in C, but this is a 1 time
    #calculation, so it does not need to be optimized).
    h_list, k_list, l_list = generate_hkl(hkl_interval, direction)
    e = generate_wt(omega_interval)
    ff_list = generate_form_factors(N_atoms_uc, atom_list, hkl_interval)


    #declare all the C variables I will need (that aren't primitives)
    cdef calg.SquareMatrix *m = calg.allocate_sq_mat(numCutOff)
    cdef np.ndarray[np.complex128_t, ndim=2] M
    cdef calg.Real **cMesh = calg.allocate_mesh(len(h_list), len(e))
    cdef calg.ComplexNum **cEvec = calg.allocate_evec(numCutOff)
    cdef calg.Real *cW = calg.allocateRealArray(numCutOff)
    cdef calg.Real *cE = dblListToCArray(e)

    #just reorginize the h,k,l points
    Q_list = np.array([h_list,k_list,l_list]).transpose()

    #benchmarking
    calc_time = time.clock()

    #Main loop - loop through all Q points
    for Qindex in range(len(h_list)):
        Q = Q_list[Qindex]

        #C implementation
        calg.createMatrix(m, clist,numCutOff, Q[0],Q[1],Q[2], magCellSize[0], magCellSize[1], magCellSize[2])

        #convert the calg.Squarematrix to a numpy array so we can find eigs
        M = sqr_mat_to_np(m)

        #find the eigs with numpy
        w, evec = np.linalg.eig(M)

        #get the form factor for this Q point
        ff = ff_list[Qindex]

        numAtoms = numCutOff

        #for testing:
        ff = 1.0

        #convert the numpy eigenvector matrix to a 2D calg.Real array in C
        convert_evec(cEvec, evec)

        #normalize the eigenvectors
        calg.normalize_evec(cEvec, clist, numCutOff)

        #convert the numpy array of eigenvalue to a C array (calg.Real)
        convert_w(cW, w)

        #Do the main calculation and fill in the appropriate row(col?) or cMesh
        calg.calc_cs(cMesh, Qindex, cEvec, cW, numCutOff, clist, cE, len(e), ff, Q[0], Q[1], Q[2], magCellSize[0], magCellSize[1], magCellSize[2], LIFETIME_VALUE)

    #convert cMesh to numpy 2D array
    mesh = np.zeros((len(h_list), len(e)))
    for i in range(len(h_list)):
        for j in range(len(e)):
            mesh[i][j] = cMesh[i][j]

    #free the memory used by variables in C
    calg.free_sq_mat(m)
    calg.free_mesh(cMesh, len(h_list), len(e))
    calg.free_evec(cEvec, numCutOff)
    calg.freeRealArray(cW)
    calg.freeRealArray(cE)
    calg.freeAtomList(clist, numCutOff)


    #calc_time = time.clock()-calc_time
    #print "generate atom list in C: ", + clist_creation, " seconds : ", 100*clist_creation/calc_time, "%"
    #print "generate various lists in python: ", gen_lists, " seconds : ", 100*gen_lists/calc_time, "%"
    #print "time creating matrices: ", create_mats, " seconds : ", 100*create_mats/calc_time, "%"
    #print "time finding eigs: ", find_eigs, " seconds : ", 100*find_eigs/calc_time, "%"
    #print "time normalizing matrices: ", norm_time, " seconds : ", 100*norm_time/calc_time, "%"
    #print "time calculating cross section points: ", cs_calc, " seconds : ", 100*cs_calc/calc_time, "%"

    print "calc time: ", time.clock()-calc_time, " seconds"
    print "Run Time: ", time.clock()-initial_time, " seconds"

    return mesh #numpy 2D array



#The following functions were copied from csection_calc, but also left in that file becuase they are useed by both
#the MyQueeny algorithm and ours.  I duplicated the code becuase the code in this file is intended to be converted to Cython
#and I assume that talking between cython and pthon is slow.




#A copy of this was also left in cesction-calc because I assume talking between python and cyton is slow
def generate_hkl(hkl_interval, direction):

    #TODO - THIS DOES NOT WORK FOR ALL CASES! only when directional components are the same or 0
    #Fixed - testing

    #What if other points are 0 (not the first or last points)?
    #For that to happen they would need to generate an odd number of points across an interval -x to x
    if hkl_interval[0] == 0:
        hkl_interval[0] = 1e-3
    elif hkl_interval[1] == 0:
        hkl_interval[1] = -1e-3

    print "hkl_interval: ", hkl_interval
    print "direction: ", direction

    #We should be multipliying these linspaces by the components, kx,ky,kz, where those are
    #normalized so that kx**2 + ky**2 + kz**2 = 1 (assuming that h**2+k**2+l**2=q**2), I'm not
    #entirely sure if this works in inverse space.
    norm = np.sqrt(direction['kx']**2 + direction['ky']**2 + direction['kz']**2)

    if direction['kx']:
        h_list = np.linspace(hkl_interval[0],hkl_interval[1],hkl_interval[2])*direction['kx']/norm
    else:
        h_list = np.zeros(hkl_interval[2])
    if direction['ky']:
        k_list = np.linspace(hkl_interval[0],hkl_interval[1],hkl_interval[2])*direction['ky']/norm
    else:
        k_list = np.zeros(hkl_interval[2])
    if direction['kz']:
        l_list = np.linspace(hkl_interval[0],hkl_interval[1],hkl_interval[2])*direction['kz']/norm
    else:
        l_list = np.zeros(hkl_interval[2])

    #print "h,k,l points: "
    #for i in range(len(h_list)):
    #    print "(",h_list[i],", ", k_list[i], ", ", l_list[i], ")"

    return h_list, k_list, l_list

#A copy of this was also left in cesction-calc because I assume talking between python and cyton is slow
def generate_wt(omega_interval):

    return np.linspace(omega_interval[0],omega_interval[1],omega_interval[2])

#A copy of this was also left in cesction-calc because I assume talking between python and cyton is slow
def generate_form_factors(N_atoms_uc, atom_list, hkl_interval):
    # Form Factor
    eval_pnts = np.linspace(hkl_interval[0],hkl_interval[1],hkl_interval[2])
    ff_list = []
    for i in range(N_atoms_uc):
        try:
            elSym = atom_list[i].label
            elSym = elements.symbol(elSym)
            elMass = atom_list[i].massNum
            if elMass == None: #The mass number will now generally be None -Tom 1/28/11
                el = elSym
            else:
                el = elSym[elMass]
            print el
            val = atom_list[i].valence
            if val != None:
                Mq = el.magnetic_ff[val].M_Q(eval_pnts)
            else:
                Mq = el.magnetic_ff[0].M_Q(eval_pnts)
        except:
            Mq = np.zeros(len(eval_pnts))
        #Test to see how things look without ff
        #Mq = np.ones(len(eval_pnts))

        ff_list = Mq
#    fake = np.ones(len(kaprange))
    return np.array(ff_list)

from csection_calc import LIFETIME_VALUE