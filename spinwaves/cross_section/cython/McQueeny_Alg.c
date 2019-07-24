#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "McQueeny_Alg.h"


//local declarations
int dbl_equal(Real d1, Real d2);

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
           //printf("freeing neighbors");
           free(list[i].neighbors);
           //printf("freeing interactions");
           free(list[i].interactions);
     }
     //printf("freeing the list");
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
Real* allocateRealArray(int size)
{
       return malloc(sizeof(Real)*size);
}
void freeRealArray(Real *list)
{
     free(list);
}

//populating the list can be done on the cython side
ComplexNum* allocateComplexArray(int size)
{
       return malloc(sizeof(ComplexNum)*size);
}
void freeComplexArray(ComplexNum *list)
{
     free(list);
}

void addToList(Atom *list, int index, int numNbrs, int *nbrs, Real *interactions, Real spin, Real pos_a, Real pos_b, Real pos_c)
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

//free the memory used by a SquareMatrix
void free_sq_mat(SquareMatrix *m)
{
     int i;
     for(i = 0; i < m->dim; i++)
     {
           free(m->mat[i]);
     }
     free(m->mat);
     free(m);
}

//allocate seperately, so we can reuse it
SquareMatrix *allocate_sq_mat(int numAtoms)
{
     int i;
     SquareMatrix *m = malloc(sizeof(SquareMatrix));
     m->mat = (ComplexNum **)malloc(sizeof(ComplexNum *)*numAtoms);
     for(i = 0; i < numAtoms; i++)
     {
           m->mat[i] = malloc(sizeof(ComplexNum)*numAtoms);
     }
     
     m->dim = numAtoms;
     return m;
}

//------------------------------Calculations------------------------------------
//This function is just converted from my Python implementation.
//numAtoms is not the size of list, but the number of atoms we want to loop
//through (the number in the magnetic cell)
SquareMatrix *createMatrix(SquareMatrix *m, Atom *list, int numAtoms, Real q0, Real q1, Real q2, int magCellSizeA, int magCellSizeB, int magCellSizeC)
{
     int i, j, lk, l;
     Real l_vec[3], l_pos[3] ,j_pos[3];
     Real tmp;
     ComplexNum cmplx_tmp;
     
     for(i = 0; i < numAtoms; i++)
     {
           for(j = 0; j < numAtoms; j++)
           {
                 //printf("assigning M[%d][%d]\n", i,j);
                 m->mat[i][j].real = 0.0;
                 m->mat[i][j].imaginary = 0.0;
                 //printf("assigned\n");
                 
                 if(i == j)
                 {
                      for(lk = 0; lk < list[i].numNeighbors; lk++)
                      {
                             //term_ij = term_ij + 2*J*spin_k
                             //This is the product of real terms and can therefore only be real
                             m->mat[i][j].real += 2*list[i].interactions[lk]*list[list[i].neighbors[lk]].spin;
                      }
                  }
         //         printf("about to do calculation for %d neighbors\n", list[i].numNeighbors);
                  //atom_j = atom_list[j]
                  j_pos[0] = list[j].pos[0];
                  j_pos[1] = list[j].pos[1];
                  j_pos[2] = list[j].pos[2];
                  for(l = 0; l < list[i].numNeighbors; l++)
                  {
           //             printf("doing calculation for neighbor %d\n", l);
                        //l_pos = list[list[i].neighbors[l]].pos;
             //           printf("lpos = (%f, %f, %f)\n", l_pos[0], l_pos[1], l_pos[2]);
             
                          //I need to use a local copy (not edit the copy in the list)
                        l_pos[0] = list[list[i].neighbors[l]].pos[0];
                        l_pos[1] = list[list[i].neighbors[l]].pos[1];
                        l_pos[2] = list[list[i].neighbors[l]].pos[2];
               //         printf("here\n");
                        
                        //l_vec = l_pos - atom_list[i].pos
                        l_vec[0] = l_pos[0] - list[i].pos[0];
                        l_vec[1] = l_pos[1] - list[i].pos[1];
                        l_vec[2] = l_pos[2] - list[i].pos[2];
                                     
                        //l_pos = l_pos/magCellSize
                        //Is this right?  I'm dividing a Real by an int, meaning
                        //I get a Real back and the fractional coordinates will
                        //be screwed up? No, I just converted to magnetic cell
                        //coordinates.
                        l_pos[0] = l_pos[0]/magCellSizeA;
                        l_pos[1] = l_pos[1]/magCellSizeB;
                        l_pos[2] = l_pos[2]/magCellSizeC;
                        
                     //   printf("here1.1");
                        
                        //l_pos = l_pos - l_pos.astype(int) #translate to first cell
                        l_pos[0] = l_pos[0] - ((int)l_pos[0]);
                        l_pos[1] = l_pos[1] - ((int)l_pos[1]);
                        l_pos[2] = l_pos[2] - ((int)l_pos[2]);
                        
                       // printf("here1.2");
                        
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
                        j_pos[0] = j_pos[0]/magCellSizeA;
                        j_pos[1] = j_pos[1]/magCellSizeB;
                        j_pos[2] = j_pos[2]/magCellSizeC;
                        
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
                        //printf("here2\n");
          //              printf("(%f,%f,%f) =? (%f,%f,%f)\n", l_pos[0], l_pos[1], l_pos[2], j_pos[0], j_pos[1], j_pos[2]);
                        if(dbl_equal(l_pos[0],j_pos[0]) && dbl_equal(l_pos[1],j_pos[1]) && dbl_equal(l_pos[2],j_pos[2]))
                        {
                             //J = (get_Jij(jnums, jmats, atom_list[i].interactions[l]))[1][1]#taking a value from diagonal matrix 
                             //J = list[i].interactions[l];
                             
                             //term_ij = term_ij - 2*sigma_j*np.sqrt(S_i*S_j)*J*np.exp(1j*np.dot(q,l_vec))
                             /*Remember:
                                        S_i = list[i].spin  (Real)
                                        S_j = list[j].spin
                                        sigma_j = S_j/abs(S_j) (+/-1.0)
                             */
                //             printf("Real equal!\n");
                             //I'll do the exponent first tmp = dot(q,l_vec)
                             tmp = q0*l_vec[0];
                             tmp += q1*l_vec[1];
                             tmp += q2*l_vec[2];
                             //tmp = e^i*tmp  e^it = cos(t) + isin(t)
                             cmplx_tmp.real = cos(tmp);                 
                             cmplx_tmp.imaginary = sin(tmp);
                             
                             //reuse tmp for the coefficient of the exponent
                             tmp = (2*list[j].spin/fabs(list[j].spin))*sqrt(list[i].spin*list[j].spin)*list[i].interactions[l];
                             
                             m->mat[i][j].real = m->mat[i][j].real - tmp*cmplx_tmp.real;
                             m->mat[i][j].imaginary = m->mat[i][j].imaginary - tmp*cmplx_tmp.imaginary;
                        }
                  }
                  //printf("scaling by 1/2\n");
                  //return Mij/2.0
                  m->mat[i][j].real = m->mat[i][j].real/2.0;
                  m->mat[i][j].imaginary = m->mat[i][j].imaginary/2.0;
                  //printf("M[%d][%d]scaled by 1/2\n", i,j);
           }
     }
     //printf("returning matrix...");
     return m;
}                                                           
                             
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


Real **allocate_mesh(int qSize, int eSize)
{
       int i;
       Real **mesh = (Real **) malloc(sizeof(Real *)*qSize);
       
       for(i = 0; i < qSize; i++)
       {
             mesh[i] = (Real *) malloc(sizeof(Real)*eSize);
       }
       
       return mesh;
}

void free_mesh(Real **mesh, int qSize, int eSize)
{
     int i;
     for(i = 0; i < qSize; i++)
     {
           free(mesh[i]);
     }
     free(mesh);
}

ComplexNum **allocate_evec(int size)
{
       int i;
       ComplexNum **evec = (ComplexNum **) malloc(sizeof(ComplexNum *) * size);
       
       for(i = 0; i < size; i++)
       {
             evec[i] = (ComplexNum *) malloc(sizeof(ComplexNum) * size);
       }
       return evec;
}

void free_evec(ComplexNum **evec, int size)
{
     int i;
     for(i = 0; i < size; i++)
     {
           free(evec[i]);
     }
     free(evec);
}
           
void normalize_evec(ComplexNum **evec, Atom *atom_list, int numAtoms)
{
     //This may not be the most efficient way of doing this, but I want to
     //copy the python code exactly for now.
     int i, n;
     Real norm;
     Real *sigma = malloc(sizeof(Real)*numAtoms);
     
     //evec2 = evec.copy()
     ComplexNum **evec2 = (ComplexNum **) malloc(sizeof(ComplexNum *)*numAtoms);
     for(i = 0; i < numAtoms; i++)
     {
           evec2[i] = (ComplexNum *) malloc(sizeof(ComplexNum)*numAtoms);
           for(n = 0; n < numAtoms; n++)
           {
                 evec2[i][n] = evec[i][n];
           }
           //for i in range(numAtoms):
            //sigma.append(s[i]/np.abs(s[i]))
           sigma[i] = atom_list[i].spin/fabs(atom_list[i].spin);
     }
     
     for(n = 0; n < numAtoms; n++)
     {
           norm = 0;
           for(i = 0; i < numAtoms; i++)
           {
                 //norm += sigma[i]*(np.abs(evec2[n][i])**2)
                 norm += sigma[i]*(evec2[n][i].real*evec2[n][i].real + evec2[n][i].imaginary*evec2[n][i].imaginary);
           }
           norm = sqrt(fabs(norm));
           
           for(i = 0; i < numAtoms; i++)
           {
                 //T[n][i] = evec2[n][i]/np.sqrt(np.abs(norm))
                 evec[n][i].real = evec2[n][i].real/norm;
                 evec[n][i].imaginary = evec2[n][i].imaginary/norm;
           }
     }
     
     //free memory
     for(i = 0; i < numAtoms; i++)
     {
           free(evec2[i]);
     }
     free(evec2);
}

/* the python code:
        evec2 = evec.copy()
        #I'm going to write this from scratch
      #  print "evec before = ", evec
        T = evec
        sigma = []
        for i in range(numAtoms):
            sigma.append(s[i]/np.abs(s[i]))
        for n in range(numAtoms):
            norm = 0
            for i in range(numAtoms):
                norm += sigma[i]*(np.abs(evec2[n][i])**2)
            for i in range(numAtoms):
      #          print "norm : ", norm
      #          print "w sign: ", np.sign(w[n])
                #T[n][i] = T[n][i]*np.sqrt((w[n]/np.abs(w[n]))/norm + 0j)
                T[n][i] = evec2[n][i]/np.sqrt(np.abs(norm))
        evec = T
*/


//evec is an array of eigenvectors, w is an array of the eigenvalues
//both should be of length numAtoms
void calc_cs(Real **mesh, int Qindex, ComplexNum **evec, Real *w, int numAtoms, Atom *atoms, Real *e, int eLen, Real ff, Real q0, Real q1, Real q2, int Na, int Nb, int Nc, Real lifeTime)
{
       int i, j;
       Real tmp, *strfac, y, lorentzian;
       //I wonder what the compiler does with constants that are too big...
       Real PI = 3.1415926535897932384626433832795028841971693993751;
       ComplexNum exp, dum;
       
       //I could probably avoid allocating this, but right now I'm just trying
       //to copy the python code as closely as possible
       strfac = (Real *) malloc(sizeof(Real)*numAtoms);
       
       for(i = 0; i < numAtoms; i++)
       {
             dum.real = 0.0;
             dum.imaginary = 0.0;
             for(j = 0; j < numAtoms; j++)
             {
                   //dum += (s[j]/np.sqrt(np.abs(s[j]))) * ff * evec[j][i]
                   tmp = atoms[j].spin/sqrt(fabs(atoms[j].spin)) * ff;
                   dum.real += tmp * evec[j][i].real;
                   dum.imaginary += tmp * evec[j][i].imaginary;
                   
                   //exp = np.exp(1.0j*(-Q[0]*basis[j][0] -Q[1]*basis[j][1] -Q[2]*basis[j][2]))
                   //for some reason these basis vectors are "translated back to the first unit cell"
                   //by subtracting the cutoffCell Dimensions.  I THINK this is becuase
                   //the actual positions I read in were from the unit cell
                   //at (Na, Nb, Nc)
                   //remember: exp(it) = cos(t) + isin(t)
                   tmp = -q0*(atoms[j].pos[0]-Na) -q1*(atoms[j].pos[1]-Nb) -q2*(atoms[j].pos[2]-Nc);
                   exp.real = cos(tmp);
                   exp.imaginary = sin(tmp);
             }
             //strfac[i] = np.abs(dum)**2
             //strfac[i] = pow(dum.real, 2) + pow(dum.imaginary, 2);
             strfac[i] = dum.real*dum.real + dum.imaginary*dum.imaginary;//quicker
             //if(Qindex == 1)
             //          printf("strfac[%d]: %f\n", i, strfac[i]);
//             printf("strfac: %f\n",strfac[i]); 
       }
       //sqw = np.zeros(len(e))
       for(i = 0; i < eLen; i++) //initialize it to zero
       {
             mesh[Qindex][i] = 0.0;
       }
       //C: sqw = mesh[Qindex]
       
       for(j = 0; j < numAtoms; j++)
       {
             for(i = 0; i < eLen; i++)
             {
                   y = lifeTime/2.0;
                   
                   //lorentzian = (1/np.pi)*y/((e[i]- w[j])**2+y**2)
                   //lorentzian = (1/PI)*y/(pow((e[i] - w[j]),2)+ pow(y,2));
                   lorentzian = (1/PI)*y/( ((e[i]-w[j])*(e[i]-w[j])) + (y*y) );//quicker
                   //if(Qindex == 1)
                   //{
                   //          printf("lorentzian[j=%d][i=%d] : %f\n", j,i,lorentzian);
                   //}
       
                   //sqw[i] = sqw[i] + strfac[j] * lorentzian 
                   mesh[Qindex][i] = mesh[Qindex][i] + strfac[j] * lorentzian;
             }
       }
//       printf("mesh[%d][10]: %f\n", Qindex, mesh[Qindex][10]);
       free(strfac);
}
                   
/*  Converted python Code:
        for i in range(numAtoms):
            dum = 0
            for j in range(numAtoms):
                dum += (s[j]/np.sqrt(np.abs(s[j]))) * ff * evec[j][i] #* np.exp(1.0j*
                        #(-Q[0]*basis[j][0] -Q[1]*basis[j][1] -Q[2]*basis[j][2]))
#                print "\n\n\nwritting to f2:\n\n\n", np.exp(1.0j*
#                        (-Q[0]*basis[j][0] -Q[1]*basis[j][1] -Q[2]*basis[j][2]))
#                print "j: ", j
#                print "basis: ", basis[j]
#                print "Q.di", (-Q[0]*basis[j][0] -Q[1]*basis[j][1] -Q[2]*basis[j][2])
                exp = np.exp(1.0j*(-Q[0]*basis[j][0] -Q[1]*basis[j][1] -Q[2]*basis[j][2]))
#                f2.write(" " + str(exp))
            strfac[i] = np.abs(dum)**2
       #     print "Q: ", Q, "    strfac[",i,"]: ", strfac[i]
#            f2.write("\n\n")
#        strfac_sums.append(0.0);
#        strfac1.append(strfac[0])
#        strfac2.append(strfac[1])
        sqw = np.zeros(len(e))
        for j in range(numAtoms):
#            strfac_sums[Qindex] +=strfac[j]
            for i in range(len(e)):
                y = LIFETIME_VALUE/2
                lorentzian = (1/np.pi)*y/((e[i]- w[j])**2+y**2)
                #test
                #strfac[j] = ff*1
                sqw[i] = sqw[i] + strfac[j] * lorentzian #he has some 'bose' factor in here
                #more test
                #sqw[i] = sqw[i] + strfac[j]
        mesh[Qindex] = sqw #* 0.5*(1 + Q[0]**2) #along x, Q_hat is just 1?
*/

int dbl_equal(Real d1, Real d2)
{
     Real min_diff = 0.000000000001;//Yes, this is arbitrary
     if((d1 - d2) < min_diff && (d2-d1) < min_diff)
            return 1;
     return 0;
}

