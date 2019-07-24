#include "complex.h"

//So that I can easilly Switch precision (mostly for the GPU, which takes a
//large performance hit with doubles vs. float(8x))
typedef float Real;

typedef struct
{
        int numNeighbors;//number of interactions (length of the following two arrays)
        int *neighbors;//an array of the indices of the other atoms that this atom interacts with
        Real *interactions;//an array of the exchange values matching ^
        Real spin; //only one dimension with this algorithm
        Real pos[3]; //atom position
} Atom;

typedef struct
{
        Real real;
        Real imaginary;
} ComplexNum;

typedef struct
{
        int dim;//this matrix is a dim*dim 2D matrix
        ComplexNum **mat; //a pointer to a 2D array of Complex numbers
} SquareMatrix;


//Allocation
Atom *createAtomList(int size);
void freeAtomList(Atom *list, int size);

//Add an atom to the atom list allocated with createAtomList(size)
void addToList(Atom *list, int index, int numNbrs, int *nbrs, Real *interactions, Real spin, Real pos_a, Real pos_b, Real pos_c);

int * allocateIntArray(int size);
void freeIntArray(int *list);

Real* allocateRealArray(int size);
void freeRealArray(Real *list);

ComplexNum* allocateComplexArray(int size);
void freeComplexArray(ComplexNum *list);

SquareMatrix *allocate_sq_mat(int numAtoms);
//fill in the matrix, m with the appropraite value for Q=(q1,q2,q3).  This is the matrix we then find the eigenvalues/vectors for
SquareMatrix *createMatrix(SquareMatrix *m, Atom *list, int numAtoms, Real q0, Real q1, Real q2, int magCellSizeA, int magCellSizeB, int magCellSizeC);
void free_sq_mat(SquareMatrix *m);

//the mesh which we would like to graph; the first index is the Q index and the second is the e index
Real **allocate_mesh(int qSize, int eSize);
void free_mesh(Real **mesh, int qSize, int eSize);

//I suppose this could have been another SquareMatrix struct
ComplexNum **allocate_evec(int size);
void free_evec(ComplexNum **evec, int size);

void normalize_evec(ComplexNum **evec, Atom *atom_list, int numAtoms);

void calc_cs(Real **mesh, int Qindex, ComplexNum **evec, Real *w, int numAtoms, Atom *atoms, Real *e, int eLen, Real ff, Real q0, Real q1, Real q2, int Na, int Nb, int Nc, Real lifeTime);

//Eigendecomposition
void eigen_decomp(SquareMatrix *m, ComplexNum *eval, SquareMatrix *evec);
void calc_eigs_on_gpu(SquareMatrix *m, ComplexNum *eval, SquareMatrix *evec, int N);



