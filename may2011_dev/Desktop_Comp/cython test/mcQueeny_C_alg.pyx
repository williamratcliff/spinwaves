import numpy as np
cimport numpy as np
import time
from time import clock
from spinwaves.spinwavecalc.readfiles import atom, readFiles
from periodictable import elements

cimport cMcQueenyAlg as calg

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
	cdef Atom *clist = calg.createAtomList(len(list))
	for i in range(len(list)):
		dblJs = []
		for j in range(list[i].interactions):
			dblJs.append(get_Jij(jnums, jmats, list[i].interactions[j]))
		calg.addTolist(clist, i, len(list[i].neighbors), intListToCArray(list[i].neighbors), dblListToCArray(dblJs), list[i].spin[spinAxis], list[i].pos[0], list[i].pos[1], list[i].pos[2])
	return clist