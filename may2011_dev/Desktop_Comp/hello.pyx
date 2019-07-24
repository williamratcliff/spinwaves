#import time
cimport numpy as np
cimport cMcQueenyAlg as calg


#Note! - I think I am still including much more in the neighbor lists than is necessary
#This is still all being sent over to the C side, so if allocation takes too long, fix that

cdef int* intListToCArray(list):
	cdef int* array = calg.allocateIntArray(len(list))
	for i in range(len(list)):
		array[i] = list[i]
	return array
	
cdef double* dblListToCArray(list):
	cdef double* array = calg.allocateDoubleList(len(list))
	for i in range(len(list)):
		array[i] = list[i]
	return array

#takes a list of python SimpleAtoms and returns a list in C	
cdef calg.Atom* createCAtomList(list, jnums, jmats, spinAxis):
	cdef calg.Atom *clist = calg.createAtomList(len(list))
	for i in range(len(list)):
		dblJs = []
		for j in range(list[i].interactions):
			dblJs.append(get_Jij(jnums, jmats, list[i].interactions[j]))
		calg.addTolist(clist, i, len(list[i].neighbors), intListToCArray(list[i].neighbors), dblListToCArray(dblJs), list[i].spin[spinAxis], 0,1,2)
	return clist

print dblListToCArray([0.234,1.234,2.534,3.45,4.53452])[0]
print dblListToCArray([0.234,1.234,2.534,3.45,4.53452])[5]
	


#tests:


#print "here"

#cdef chello.Atom a
#chello.populateAtom(&a, 2)


#print "here2"
#cdef chello.SquareMatrix m = chello.getD2Matrix()
#print "here3"

#print "dim = "
#print (m).dim

#cdef myFunc(np.ndarray[np.complex128_t, ndim=2] Mij):
#	Mij[0][0] = Mij[1][1]