#include <stdio.h>
#include <stdlib.h>
#include "McQueeny_Alg.h"
#include <math.h>


//local declarations
int dbl_equal(double d1, double d2);

//-------------------------Creating the list on the C side----------------------

Atom *createAtomList(int size)
{
     return (Atom *) malloc(sizeof(Atom)*size);
}

void freeAtomList(Atom *list, int size)
{
     int i;
     for(i = 0; i < size; i++)
     {
           free(&list[i].neighbors);
           free(&list[i].interactions);
     }
     free(list);
}

//populating the list can be done on the cython side
int * allocateIntArray(int size)
{
     return malloc(size*sizeof(int));
}
void freeIntArray(int *list)
{
     free(list);
}

//populating the list can be done on the cython side
double* allocateDoubleArray(int size)
{
       return malloc(sizeof(double)*size);
}
void freeDoubleArray(double *list)
{
     free(list);
}

void addToList(Atom *list, int index, int numNbrs, int *nbrs, double *interactions, double spin, double pos_a, double pos_b, double pos_c)
{
     list[index].numNeighbors = numNbrs;
     list[index].neighbors = nbrs;//an array of the indices of the other atoms that this atom interacts with
     list[index].interactions = interactions;//an array of the exchange values matching ^
     list[index].spin = spin;
     //I could do this with a pointer, but this way I minimize memory allocations
     list[index].pos[0] = pos_a;
     list[index].pos[1] = pos_b;
     list[index].pos[2] = pos_c;
     
}

//------------------------------Calculations------------------------------------
//This function is just converted from my Python implementation.
//numAtoms is not the size of list, but the number of atoms we want to loop
//through (the number in the magnetic cell)
SquareMatrix *createMatrix(Atom *list, int numAtoms, double q1, double q2, double q3, int magCellA, int magCellB, int magCellC)
{
     int i, j, lk, l;
     double l_vec[3], *l_pos ,j_pos[3], poss_diff[3];
     double tmp;
     ComplexNum cmplx_tmp;
     //ComplexNum term_ij;
     SquareMatrix *m = malloc(sizeof(SquareMatrix));
     m->mat = malloc(sizeof(ComplexNum)*numAtoms*numAtoms);
     m->dim = numAtoms;
     
     for(i = 0; i < numAtoms; i++)
     {
           for(j = 0; j < numAtoms; j++)
           {
                 m->mat[i][j].real = 0.0;
                 m->mat[i][j].imaginary = 0.0;
                 
                 if(i == j)
                 {
                      for(lk = 0; lk < list[i].numNeighbors; lk++)
                      {
                             //term_ij = term_ij + 2*J*spin_k
                             //This is the product of real terms and can therefore only be real
                             m->mat[i][j].real += 2*list[i].interactions[lk]*list[list[i].neighbors[lk]].spin;
                      }
                  }
                  for(l = 0; l < list[i].numNeighbors; l++)
                  {
                        l_pos = list[list[i].neighbors[l]].pos;
                        //l_pos[0] = list[list[i].neighbors[l]].pos[0];
                        //l_pos[1] = list[list[i].neighbors[l]].pos[1];
                        //l_pos[2] = list[list[i].neighbors[l]].pos[2];
                        
                        //l_vec = l_pos - atom_list[i].pos
                        l_vec[0] = l_pos[0] - list[i].pos[0];
                        l_vec[1] = l_pos[1] - list[i].pos[1];
                        l_vec[2] = l_pos[2] - list[i].pos[2];
                        
                        //l_pos = l_pos/magCellSize
                        l_pos[0] = l_pos[0]/magCellA;
                        l_pos[1] = l_pos[1]/magCellB;
                        l_pos[2] = l_pos[2]/magCellC;
                        
                        //l_pos = l_pos - l_pos.astype(int) #translate to first cell
                        l_pos[0] = l_pos[0] - ((int)l_pos[0]);
                        l_pos[1] = l_pos[1] - ((int)l_pos[1]);
                        l_pos[2] = l_pos[2] - ((int)l_pos[2]);
                        
                        //#in case l_pos is negative:
                        //l_pos = np.ceil(-l_pos) + l_pos
                        //is this really necessary?  I don't see how l_pos could be negative.
                        //l_vec obviously could be, but l_pos is just a position.
                        //I feel like I put that line in the python code for 
                        //a reason though, so I'll leave it here.
                        //It shouldn't be too expensive
                        l_pos[0] = ceil(-l_pos[0]) + l_pos[0];
                        l_pos[1] = ceil(-l_pos[1]) + l_pos[1];
                        l_pos[2] = ceil(-l_pos[2]) + l_pos[2];
                        
                        //repeat this for the j_pos
                        //j_pos = atom_j.pos/magCellSize
                        j_pos[0] = j_pos[0]/magCellA;
                        j_pos[1] = j_pos[1]/magCellB;
                        j_pos[2] = j_pos[2]/magCellC;
                        
                        //j_pos = j_pos - j_pos.astype(int)
                        j_pos[0] = j_pos[0] - ((int)j_pos[0]);
                        j_pos[1] = j_pos[1] - ((int)j_pos[1]);
                        j_pos[2] = j_pos[2] - ((int)j_pos[2]);
                        
                        //j_pos = np.ceil(-j_pos) + j_pos
                        j_pos[0] = ceil(-j_pos[0]) + j_pos[0];
                        j_pos[1] = ceil(-j_pos[1]) + j_pos[1];
                        j_pos[2] = ceil(-j_pos[2]) + j_pos[2];
                        
                        //pos_diff = (l_pos)-(j_pos)#converted to magcell pos
                        //poss_diff[0] = l_pos[0] - j_pos[0];
                        //poss_diff[1] = l_pos[1] - j_pos[1];
                        //poss_diff[2] = l_pos[2] - j_pos[2];
                        
                        //if ((np.abs(pos_diff)<min_diff).all()):
                        //if(abs(pos_diff)
                        if(dbl_equal(l_pos[0],j_pos[0]) && dbl_equal(l_pos[1],j_pos[1]) && dbl_equal(l_pos[2],j_pos[2]))
                        {
                             //J = (get_Jij(jnums, jmats, atom_list[i].interactions[l]))[1][1]#taking a value from diagonal matrix 
                             //J = list[i].interactions[l];
                             
                             //term_ij = term_ij - 2*sigma_j*np.sqrt(S_i*S_j)*J*np.exp(1j*np.dot(q,l_vec))
                             /*Remember:
                                        S_i = list[i].spin  (double)
                                        S_j = list[j].spin
                                        sigma_j = S_j/abs(S_j) (+/-1.0)
                             */
                             //I'll do the exponent first tmp = dot(q,l_vec)
                             tmp = q1*l_vec[0];
                             tmp += q2*l_vec[1];
                             tmp += q3*l_vec[2];
                             //tmp = e^i*tmp  e^it = cos(t) + isin(t)
                             cmplx_tmp.real = cos(tmp);                 
                             cmplx_tmp.imaginary = sin(tmp);
                             
                             //reuse tmp for the coefficient of the exponent
                             tmp = (2*list[j].spin/abs(list[j].spin))*sqrt(list[i].spin*list[j].spin)*list[i].interactions[l];
                             
                             m->mat[i][j].real = m->mat[i][j].real - tmp*cmplx_tmp.real;
                             m->mat[i][j].imaginary = m->mat[i][j].imaginary - tmp*cmplx_tmp.imaginary;
                        }
                  }
                  //return Mij/2.0
                  m->mat[i][j].real = m->mat[i][j].real/2.0;
                  m->mat[i][j].imaginary = m->mat[i][j].imaginary/2.0;
           }
     }
     
     return m;
                                                             
                             
/*  the converted python code:
     
for i in range(natoms):
        for j in range(natoms):
            term_ij = 0.0 + 0j
            #print "i: ", i, "  j: ", j
            #print "atom_i pos: ", atom_list[i].pos
            if i == j:#delta_ij term (diaginal elements only)
                for lk in range(len(atom_list[i].neighbors)):
                    #assuming Jij matrix is diagonal, take one of the diagonal terms:
                    J = (get_Jij(jnums, jmats,atom_list[i].interactions[lk]))[1][1]
                    spin_k = atom_list[atom_list[i].neighbors[lk]].spin[spin_comp]
                    term_ij = term_ij + 2*J*spin_k#*sigma_k*S_k
            S_i = np.abs(atom_list[i].spin[spin_comp])
            spin_j = atom_list[j].spin[spin_comp]
            sigma_j = spin_j/np.abs(spin_j)
            S_j = np.abs(spin_j)
            #now for each atom correspondingto atom j in cell l sum J * exp(i*Q.l)
            #I'll loop through all the interactions, but only pick out those that bond
            #corresponding atoms in different(or the same) cells, ie. their positions - any integer
            #component are the same.
            atom_j = atom_list[j]
            for l in range(len(atom_list[i].neighbors)):
                atom_l = atom_list[atom_list[i].neighbors[l]]
                #Use unit cell lengtsh to construct l vecter
                l_pos = atom_l.pos
                #print "l_pos: ", l_pos
                #First we assumed it had to be a vector with integer components:
                #l_vec = l_pos.astype(int)-(atom_list[i].pos).astype(int)#firstCellPos
                #However, that ^ often produced imaginary numbers in the matrix
                l_vec = l_pos - atom_list[i].pos
                #use magnetic cell lengths to determine if atom_l is has the same position
                l_pos = l_pos/magCellSize
                l_pos = l_pos - l_pos.astype(int) #translate to first cell
                #in case l_pos is negative:
                l_pos = np.ceil(-l_pos) + l_pos
                #same for j_pos
                j_pos = atom_j.pos/magCellSize
                j_pos = j_pos - j_pos.astype(int)
                j_pos = np.ceil(-j_pos) + j_pos
                
                pos_diff = (l_pos)-(j_pos)#converted to magcell pos
                if ((np.abs(pos_diff)<min_diff).all()):
                    #print "l_vec: ", l_vec
                    J = (get_Jij(jnums, jmats, atom_list[i].interactions[l]))[1][1]#taking a value from diagonal matrix  <---------Here
                    #print "sinusoidal comp: ", 2*sigma_j*np.sqrt(S_i*S_j)*J*np.exp(1j*np.dot(q,l_vec))
                    term_ij = term_ij - 2*sigma_j*np.sqrt(S_i*S_j)*J*np.exp(1j*np.dot(q,l_vec))
                    #print "term[",i,"][",j,"]: ", term_ij
   #             elif (i != j):
                    #for testing:
    #                print ""
                    #print "\n\n\nOMG!!!  I shouldn't be here for the simple ferro/antiferro cases!\n\n\n"
            #print "term_ij: ", term_ij
            Mij[i][j] = term_ij
            #print "Mij: ", Mij
    #print "Q: ", q
    #print "M:\n", Mij
    #this creates eigenvalues that are different from ours by a factor of 2.  For now I'll
    #I'll just divide that out.
    return Mij/2.0
    */
}


int dbl_equal(double d1, double d2)
{
     double min_diff = 0.0000000000001;
     if((d1 - d2) < min_diff && (d2-d1) < min_diff)
            return 1;
     return 0;
}


//void freeMatrix(SquareMatrix *m)//frees the memory this matrix is using
//{
 //    free(m->mat);
  //   free(m);
//}


//----------------------------Some Cython Tests---------------------------------

//A test
void populateAtom(Atom *a, int numNbrs)
{
     a->numNeighbors = numNbrs;
}

// A test
SquareMatrix getD2Matrix()
{
     SquareMatrix *s = malloc(sizeof(SquareMatrix));
     s->dim = 2;
     s->mat = (ComplexNum**)malloc(sizeof(ComplexNum[2][2]));
     s->mat[0][0].real = 1;
     s->mat[0][1].real = 2;
     s->mat[1][0].real = 3;
     s->mat[1][1].real = 4;
     return *s;
}

void printList(int *list, int size)
{
     int i;
     for(i = 0; i < size; i++)
     {
           printf("%d) %d\n", i, list[i]);
     }
}
     

void doNothing(int *list)
{
     printf("nothing");
}

