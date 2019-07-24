#include "complex.h"

typedef struct
{
        int numNeighbors;//number of interactions (length of the following two arrays)
        int *neighbors;//an array of the indices of the other atoms that this atom interacts with
        double *interactions;//an array of the exchange values matching ^
        double spin; //only one dimension with this algorithm
        double pos[3]; //atom position
} Atom;

typedef struct
{
        double real;
        double imaginary;
} ComplexNum;

typedef struct
{
        int dim;
        ComplexNum **mat; //a pointer to a 2D array of Complex numbers
} SquareMatrix;


//Allocation
Atom *createAtomList(int size);
void freeAtomList(Atom *list, int size);

void addToList(Atom *list, int index, int numNbrs, int *nbrs, double *interactions, double spin, double pos_a, double pos_b, double pos_c);

int * allocateIntArray(int size);
void freeIntArray(int *list);

double* allocateDoubleArray(int size);
void freeDoubleArray(double *list);


//test
void printList(int *list, int size);
void doNothing(int* list);


