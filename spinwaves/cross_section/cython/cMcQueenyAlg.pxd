cdef extern from "McQueeny_Alg.h":
	
	ctypedef double Real
	
	ctypedef struct Atom:
		pass
		#int numNeighbors#number of interactions (length of the following two arrays)
		#int *neighbors#and array of the indices of the other atoms that this atom interacts with
		#float *interactions#an array of the exchange values matching ^
		#float spin[3]

	ctypedef struct ComplexNum:
		Real real
		Real imaginary
	ctypedef struct SquareMatrix:
		int dim
		ComplexNum **mat #a 2D array of Complex numbers

	#//Allocation
	Atom *createAtomList(int size)
	void freeAtomList(Atom *list, int size)
	
	void addToList(Atom *list, int index, int numNbrs, int *nbrs, Real *interactions, Real spin, Real pos_a, Real pos_b, Real pos_c)
	
	int * allocateIntArray(int size)
	void freeIntArray(int *list)
	
	Real* allocateRealArray(int size)
	void freeRealArray(Real *list)

	ComplexNum* allocateComplexArray(int size)
	void freeComplexArray(ComplexNum *list)

	SquareMatrix *allocate_sq_mat(int numAtoms)
	SquareMatrix *createMatrix(SquareMatrix *m, Atom *list, int numAtoms, Real q0, Real q1, Real q2, int magCellSizeA, int magCellSizeB, int magCellSizeC)
	void free_sq_mat(SquareMatrix *m)

	Real **allocate_mesh(int qSize, int eSize)
	void free_mesh(Real **mesh, int qSize, int eSize)

	ComplexNum **allocate_evec(int size)
	void free_evec(ComplexNum **evec, int size)

	void normalize_evec(ComplexNum **evec, Atom *atom_list, int numAtoms)

	void calc_cs(Real **mesh, int Qindex, ComplexNum **evec, Real *w, int numAtoms, Atom *atoms, Real *e, int eLen, Real ff, Real q0, Real q1, Real q2, int Na, int Nb, int Nc, Real lifeTime)
	
	
	#//test
	void printList(int *list, int size)
	void doNothing(int* list)