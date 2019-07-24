cdef extern from "McQueeny_Alg.h":
	ctypedef struct Atom:
		pass
		#int numNeighbors#number of interactions (length of the following two arrays)
		#int *neighbors#and array of the indices of the other atoms that this atom interacts with
		#float *interactions#an array of the exchange values matching ^
		#float spin[3]

	ctypedef struct ComplexNum:
		double real
		double imaginary
	ctypedef struct SquareMatrix:
		int dim
		ComplexNum **mat #a 2D array of Complex numbers

	#//Allocation
	Atom *createAtomList(int size)
	void freeAtomList(Atom *list, int size)
	
	void addToList(Atom *list, int position, int numNbrs, int *nbrs, double *interactions, double spin, double pos_a, double pos_b, double pos_c)
	
	int * allocateIntArray(int size)
	void freeIntList(int *list)
	
	double* allocateDoubleArray(int size)
	void freeDoubleList(double *list)
	
	
	#//test
	void printList(int *list, int size)
	void doNothing(int* list)