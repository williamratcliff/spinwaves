#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "McQueeny_Alg.h"


//local declarations
int dbl_equal(Real d1, Real d2);
void tred2(Real *V, Real *d, Real *e, int n);
void tql2(Real *V, Real *d, Real *e, int n);

//Putting it togethor
//really, the e's could be computed in here, along withthe (h,k,l)s
void calc_csection(Real *h, Real *k, Real *l, int numPoints, Atom *atoms, int numAtoms, Real *ff_list, Real *e, int e_len, int magCellSizeA, int magCellSizeB, int magCellSizeC, Real** mesh, Real LifeTime)
{
     int i,j;
     SquareMatrix *m = allocate_sq_mat(numAtoms);
     SquareMatrix *evec = allocate_sq_mat(numAtoms);
     ComplexNum *eval = allocateComplexArray(numAtoms);
     Real *w = allocateRealArray(numAtoms);
     
     for(i = 0; i < numPoints; i++)
     {
           createMatrix(m, atoms, numAtoms, h[i], k[i], l[i], magCellSizeA, magCellSizeB, magCellSizeC);
           eigen_decomp(m, eval, evec);
           normalize_evec(evec->mat, atoms, numAtoms);
           //the evals should be real
           for(j = 0; j < numAtoms; j++)
           {
                 w[j] = eval[j].real;
           }
           calc_cs(mesh, i, evec->mat, w, numAtoms, atoms, e, e_len, ff_list[i], h[i], k[i], l[i], magCellSizeA, magCellSizeB, magCellSizeC, LifeTime);
     }
     freeComplexArray(eval);
     free_sq_mat(m);
     free_sq_mat(evec);
     freeRealArray(w);
}



//-----------------------Eigenvalue/vector Code---------------------------------

/*I am taking the releveant code from TNT and adapting it to our problem in C.
I am taking only the code necessary to find the eigenvectors of a hermitian
matrix.  The TNT code is only for real matrices, but luckilly a Hermitian
matrix can be converted to a real matrix to find the eigenvalues/vectors.
*/


/*
This function will find the eigenvalues of matrix m and put them into the
array eval.  The eigenvectors will be put into evec.

-m must be Hermitian.
*/
void eigen_decomp(SquareMatrix *m, ComplexNum *eval, SquareMatrix *evec)
{
     int i,j;
     int N = m->dim;
     int n = N*2;//The TNT code uses this variable name
     
     //first I need to construct the larger, real matrix.
     
     /* From Numerical Recipes in C, 11.4:
             
             Our complex matrix, C, can be represented like this:
             
                 C = A + iB
             
             Our Eigenvalue problem then becomes:
                 
                 (A + iB)*(u + iv) = l(u + iv)
                 
                 | A  -B | | u |     | u |
                 |       | |   | = l |   |
                 | B   A | | v |     | v |
                 
             The eigenvalues in l will then come in pairs: l1,l1,l2,l2,l3,l3...
             and the eigenvectors of C will be u +iv and i(u+iv)
     */
     
     //I am allocating a single, linear array for this 2D array
     //To index it:  A[i][j] = A[i*N+j], where A is an N*N matrix
     //I'll keep the same names as the TNT code.  They use V for the
     //eigenvectors, but initially copy the matrix A to V, and then don't touch
     //A.  
     Real *V = (Real *) malloc(sizeof(Real)*n*n);
     Real *d = (Real *) malloc(sizeof(Real)*n);// real part of eigenvalues
     Real *e = (Real *) malloc(sizeof(Real)*n);// imag part of eigenvalues
     
     //populate A
     for(i = 0; i < N; i++)
     {
           for(j = 0; j < N; j++)
           {
                 //printf("putting %f at %d\n", (m->mat)[i][j].real, i*n + j);
                 V[i*n + j] = (m->mat)[i][j].real;
                 V[i*n + (N+j)] = -(m->mat)[i][j].imaginary;
                 V[(N+i)*n + j] = (m->mat)[i][j].imaginary;
                 V[(N+i)*n + (N+j)] = (m->mat)[i][j].real;
           }
     }
     
     //test
     //n = 3;
     //V = (Real *) malloc(sizeof(Real)*n*n);
     //V[0] = 1.0;
     //V[1] = 2.0;
     //V[2] = 3.0;
     //V[3] = 2.0;
     //V[4] = 1.5;
     //V[5] = -2.5;
     //V[6] = 3.0;
     //V[7] = -2.5;
     //V[8] = 0.5;
     
//     printf("The imag Matrix:\n");
//     for(i = 0; i < N; i++)
//     {
//           for(j = 0; j < N; j++)
//           {
//                 printf("%f + %fi     ", m->mat[i][j].real, m->mat[i][j].imaginary);
//           }
//           printf("\n");
//     }
     
     //hardcode my real matrix here for testing
    // n = 4;
//     V = (Real *) malloc(sizeof(Real)*n*n);
//     V[0] = 1;
//     V[1] = 5;
//     V[2] = 3;
//     V[3] = 7;
//     V[4] = 5;
//     V[5] = 4;
//     V[6] = 8;
//     V[7] = 3;
//     V[8] = 3;
//     V[9] = 8;
//     V[10] = 1;
//     V[11] = 7;
//     V[12] = 7;
//     V[13] = 3;
//     V[14] = 7;
//     V[15] = 5;
     
     
//     printf("The real Matrix:\n");
//     for(i = 0; i < n; i++)
//     {
//           for(j = 0; j < n; j++)
//           {
//                 printf("%f  ", V[i*n + j]);
//           }
//           printf("\n");
//     }
     
     //Tridiagonalize
     tred2(V, d, e, n);
     
     
//     printf("\n\n\n---------------------------- Diagonalizing ------------------------------\n\n");
     
     //Diagonalize
     tql2(V, d, e, n);
     
     //Now reconstruct a complex answer
     //numpy can be used to figure out this ordering, but right now the eigs for
     //my big matrix are wrong.
     
     //test
//     printf("\n\neigs:\n");
//     for(i = 0; i < n; i++)
//     {
//           //printf("eval: %f + %fi\n", d[i], e[i]);
//           //printf("evec row:\n");
//           for(j = 0; j < n; j++)
//           {
//                 printf("%f ", V[i*n + j]);
//           }
//           printf("\n");
//     }
//     for(i = 0; i < n; i++)
//     {
//           printf("eval: %f + %fi\n", d[i], e[i]);
//     }
     
     //Now to construct the complex eigenvectors
     
     //The vectors/values are sorted by the eigenvalue, so redundant pairs will
     //appear next to each other.  Therefore, I just take every other
     //vector/value
     //construct the complex eigenvector matric and value array
     for(j = 0; j < n; j+=2)//column/value
     {
             eval[j/2].real = d[j];
             eval[j/2].imaginary = e[j];
             
             for(i = 0; i < N; i++)//row
             {
                   //printf("mat[%d][%d] = %f + %fi\n", i, j, V[i*n + j], V[(i+N)*n+j]);
                   evec->mat[i][j/2].real = V[i*n + j];
                   evec->mat[i][j/2].imaginary = V[(i+N)*n + j];
             }
     }
     free(V);
     free(d);
     free(e);
     
     //print the complex evecs for testing
//     printf("\n\ncomplex eigenvectors:\n");
//     for(i = 0; i < N; i++)
//     {
//           for(j = 0; j < N; j++)
//           {
//                 printf("%f +%fi   ", evec->mat[i][j].real, evec->mat[i][j].imaginary);
//           }
//           printf("\n");
//     }
}



//Adapted from the TNT code (Eigenvalue class) in jama_eig.h
void tred2(Real *V, Real *d, Real *e, int n)
{
//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

 //   int tmpi, tmpj;//for testing
    int i,j,k;
    Real scale, h, f, g, hh;
    
    for (j = 0; j < n; j++)
    {
        //d[j] = V[n-1][j];
        d[j] = V[(n-1)*n+j];
    }
    
//    printf("e initial: %E %E %E %E\n",e[0],e[1],e[2],e[3]); 

    // Householder reduction to tridiagonal form.

    for (i = n-1; i > 0; i--)
    {
//    printf("\nIn for loop(top):\n");
//     for(tmpi = 0; tmpi < n; tmpi++)
//     {
//           for(tmpj = 0; tmpj < n; tmpj++)
//           {
//                 printf("%f  ", V[tmpi*n + tmpj]);
 //          }
  //         printf("\n");
 //    }

        // Scale to avoid under/overflow.
        
        //Real scale = 0.0;
        //Real h = 0.0;
        
        //I haven't looked through all the code for all the uses of these, but
        //I think it's very unlikely that they can't be overwritten for each
        //iteration.
        scale = 0.0;
        h = 0.0;
        
//        printf("d before scale calc: %E %E %E %E\n",d[0],d[1],d[2],d[3]); 
        for (k = 0; k < i; k++)
        {
            scale = scale + fabs(d[k]);
            //printf("d[%d] = %f   | abs = %f\n", k, d[k], abs(d[k]));
            //printf("Scale = %f\n", scale);
        }
        if (scale == 0.0)
        {
                  
//        printf("\nscale =0:\n");
//         for(tmpi = 0; tmpi < n; tmpi++)
 //        {
 //              for(tmpj = 0; tmpj < n; tmpj++)
 //              {
 //                    printf("%f  ", V[tmpi*n + tmpj]);
 //              }
 //              printf("\n");
 //        }
  //                
//           printf("scale 0:\n");
//           t = i;
//           for(i = 0; i < n; i++)
//           {
//            for(j = 0; j < n; j++)
//            {
//                  printf("%f  ", V[i*n + j]);
//            }
//            printf("\n");
//           }
//           i = t;
           
           
            e[i] = d[i-1];
//            printf("e after first assign: %E %E %E %E\n",e[0],e[1],e[2],e[3]); 
            for (j = 0; j < i; j++)
            {
                //d[j] = V[i-1][j];
                d[j] = V[(i-1)*n +j];
                
                //V[i][j] = 0.0;
                V[i*n +j] = 0.0;
                
                //V[j][i] = 0.0;
                V[j*n +i] = 0.0;
            }
        } else 
        {
               
//        printf("\nscale not 0:\n");
//        printf("scale = %E\n", scale);
//         for(tmpi = 0; tmpi < n; tmpi++)
//         {
 //              for(tmpj = 0; tmpj < n; tmpj++)
  //             {
   //                  printf("%f  ", V[tmpi*n + tmpj]);
   //            }
  //             printf("\n");
  //       }
  //       printf("\nf: %f    g: %f\n", f, g);
  //       printf("d0.5: ");
  //       for(tmpi = 0; tmpi < n; tmpi++)
  //       {
  //                printf("%f  ", d[tmpi]);
  //       }
  //       printf("\n");
        
            // Generate Householder vector.
            for (k = 0; k < i; k++)
            {
//                printf("d before scale: %E %E %E %E\n",d[0],d[1],d[2],d[3]); 
                d[k] /= scale;
//                printf("d after scale: %E %E %E %E\n",d[0],d[1],d[2],d[3]);
                h += d[k] * d[k];
            }

            //Real f = d[i-1];
            //Real g = sqrt(h);
            f = d[i-1];
            g = sqrt(h);

//         printf("d2: ");
 //        for(tmpi = 0; tmpi < n; tmpi++)
 //        {
 //                 printf("%f  ", d[tmpi]);
 //        }
 //        printf("\n");
//         printf("\nf: %f    g: %f", f, g);

            if (f > 0)
            {
                  g = -g;
            }
            
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            
            
            
            for (j = 0; j < i; j++)
            {
                e[j] = 0.0;
            }
//            printf("e : %E %E %E %E\n",e[0],e[1],e[2],e[3]); 
//            printf("d: ");
//         for(tmpi = 0; tmpi < n; tmpi++)
//         {
//                  printf("%f  ", d[tmpi]);
//         }
//         printf("\n");
            // Apply similarity transformation to remaining columns.

            for (j = 0; j < i; j++)
            {
                f = d[j];
                //V[j][i] = f;
                V[j*n+i] = f;
                //g = e[j] + V[j][j] * f;
                g = e[j] + V[j*n + j] * f;
                
                for (k = j+1; k <= i-1; k++)
                {
                    //g += V[k][j] * d[k];
                    g += V[k*n+j] * d[k];
                    
                    //e[k] += V[k][j] * f;
                    e[k] += V[k*n+j] * f;
                }
                e[j] = g;
            }
//            printf("e : %E %E %E %E\n",e[0],e[1],e[2],e[3]); 
//            printf("here1\n");
//         for(tmpi = 0; tmpi < n; tmpi++)
//         {
//               for(tmpj = 0; tmpj < n; tmpj++)
//               {
//                     printf("%f  ", V[tmpi*n + tmpj]);
//               }
//               printf("\n");
//         }
            
            
            
            f = 0.0;
            
            for (j = 0; j < i; j++)
            {
                e[j] /= h;
                f += e[j] * d[j];
            }
//            printf("e2.25: %E %E %E %E\n",e[0],e[1],e[2],e[3]);
            //Real hh = f / (h + h);
            hh = f / (h + h);
//            printf("\nhh: %E\n", hh);
            
            for (j = 0; j < i; j++)
            {
                e[j] -= hh * d[j];
            }
//            printf("d2.5: %E %E %E %E\n",d[0],d[1],d[2],d[3]);
//           printf("e2.5: %E %E %E %E\n",e[0],e[1],e[2],e[3]);
            for (j = 0; j < i; j++)
            {
                f = d[j];
                g = e[j];
                
                for (k = j; k <= i-1; k++) 
                {
                    //V[k][j] -= (f * e[k] + g * d[k]);
                    V[k*n+j] -= (f * e[k] + g * d[k]);
                }
                
                //d[j] = V[i-1][j];
                d[j] = V[(i-1)*n+j];
                
                //V[i][j] = 0.0;
                V[i*n+j] = 0.0;
            }
//            printf("d3: %E %E %E %E\n",d[0],d[1],d[2],d[3]);
        }
        d[i] = h;
    
//    printf("\nIn for loop(bottom):\n");
//     for(tmpi = 0; tmpi < n; tmpi++)
//     {
//           for(tmpj = 0; tmpj < n; tmpj++)
//           {
//                 printf("%f  ", V[tmpi*n + tmpj]);
//           }
//           printf("\n");
//     }
//      printf("d at end of loop: %E %E %E %E\n",d[0],d[1],d[2],d[3]);
    
    }
    
    
//    printf("Before Accumulating Transforms:\n");
//     for(i = 0; i < n; i++)
//     {
//           for(j = 0; j < n; j++)
 //          {
  //               printf("%f  ", V[i*n + j]);
 //          }
 //          printf("\n");
 //    }
    // Accumulate transformations.
    
    for (i = 0; i < n-1; i++)
    {
        //V[n-1][i] = V[i][i];
        V[(n-1)*n +i] = V[i*n +i];
        
        //V[i][i] = 1.0;
        V[i*n +i] = 1.0;
        
        //Real h = d[i+1];
        h = d[i+1];
        
        if (h != 0.0)
        {
              for (k = 0; k <= i; k++)
              {
                  //d[k] = V[k][i+1] / h;
                  d[k] = V[k*n +i+1] / h;
              }
              
              for (j = 0; j <= i; j++)
              {
                  //Real g = 0.0;
                  g = 0.0;
                  
                  for (k = 0; k <= i; k++)
                  {
                      //g += V[k][i+1] * V[k][j];
                      g += V[k*n +i+1] * V[k*n +j];
                  }
                  for (k = 0; k <= i; k++)
                  {
                      //V[k][j] -= g * d[k];
                      V[k*n+j] -= g * d[k];
                  }
              }
        }
        
        for (k = 0; k <= i; k++)
        {
            //V[k][i+1] = 0.0;
            V[k*n + i+1] = 0.0;
        }
    }
    
 //   printf("After Accumulating Transforms:\n");
//     for(i = 0; i < n; i++)
//     {
//           for(j = 0; j < n; j++)
//           {
//                 printf("%f  ", V[i*n + j]);
//           }
//           printf("\n");
//     }
    
    for (j = 0; j < n; j++)
    {
        //d[j] = V[n-1][j];
        d[j] = V[(n-1)*n +j];
        
        //V[n-1][j] = 0.0;
        V[(n-1)*n +j] = 0.0;
    }
    
    //V[n-1][n-1] = 1.0;
    V[(n-1)*n +n-1] = 1.0;
    e[0] = 0.0;
    
    
//     printf("\n\nThe tridiagonal Matrix:\n");
//     for(i = 0; i < n; i++)
//     {
//           for(j = 0; j < n; j++)
//           {
//                 printf("%f  ", V[i*n + j]);
//           }
//           printf("\n");
//     }
    
}




//Adapted from function in Eigenvalue class in TNT code

// Symmetric tridiagonal QL algorithm.
   
void tql2(Real *V, Real *d, Real *e, int n)
{
     //  This is derived from the Algol procedures tql2, by
     //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
     //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
     //  Fortran subroutine in EISPACK.
     
     
//     int tmpi, tmpj; //testing
     int i, l, m, iter, k, j;
     Real f, tst1, eps, g, p, r, dl1, h, c, c2, c3, el1, s, s2;
     
//     printf("d1: %E %E %E %E\n",d[0],d[1],d[2],d[3]);
//     printf("e1: %E %E %E %E\n",e[0],e[1],e[2],e[3]);
//     printf("V1:\n");
//         for(tmpi = 0; tmpi < n; tmpi++)
//     {
//           for(tmpj = 0; tmpj < n; tmpj++)
//           {
//                 printf("%f  ", V[tmpi*n + tmpj]);
//           }
//           printf("\n");
//     }
//     printf("\n\n");
     
     for (i = 1; i < n; i++)
     {
         e[i-1] = e[i];
     }
     
     e[n-1] = 0.0;
     
     //Real f = 0.0;
     //Real tst1 = 0.0;
     //Real eps = pow(2.0,-52.0);
     f = 0.0;
     tst1 = 0.0;
     eps = pow(2.0, -52.0);
     
//     printf("d2: %E %E %E %E\n",d[0],d[1],d[2],d[3]);
//     printf("e2: %E %E %E %E\n",e[0],e[1],e[2],e[3]);
//     printf("V2:\n");
 //         for(tmpi = 0; tmpi < n; tmpi++)
 //    {
 //          for(tmpj = 0; tmpj < n; tmpj++)
 //          {
 //                printf("%f  ", V[tmpi*n + tmpj]);
 //          }
 //          printf("\n");
 //    }
 //    printf("\n\n");
     
     
     for (l = 0; l < n; l++)
     {
         // Find small subdiagonal element
         
         tst1 = fmax(tst1,fabs(d[l]) + fabs(e[l]));
         //int m = l;
         m = l;
         
         // Original while-loop from Java code
         
         while (m < n)
         {
              if (fabs(e[m]) <= eps*tst1)
              {
                   break;
              }
              m++;
         }
         
         // If m == l, d[l] is an eigenvalue,
         // otherwise, iterate.
         
         if (m > l)
         {
               //int iter = 0;
               iter = 0;
               
               do 
               {
                   iter = iter + 1;  // (Could check iteration count here.)
                   
                   // Compute implicit shift
                   
                   //Real g = d[l];
                   //Real p = (d[l+1] - g) / (2.0 * e[l]);
                   //Real r = hypot(p,1.0);
                   g = d[l];
                   p = (d[l+1] - g) / (2.0 * e[l]);
                   r = hypot(p,1.0);
                   
                   
                   if (p < 0)
                   {
                         r = -r;
                   }
                   
                   d[l] = e[l] / (p + r);
                   d[l+1] = e[l] * (p + r);
                   
                   //Real dl1 = d[l+1];
                   //Real h = g - d[l];
                   dl1 = d[l+1];
                   h = g - d[l];
                   
                   for (i = l+2; i < n; i++)
                   {
                       d[i] -= h;
                   }
                   
                   f = f + h;
                   
 //                  printf("d3: %E %E %E %E\n",d[0],d[1],d[2],d[3]);
 //    printf("e3: %E %E %E %E\n",e[0],e[1],e[2],e[3]);
 //    printf("V3:\n");
 //         for(tmpi = 0; tmpi < n; tmpi++)
 //    {
 //          for(tmpj = 0; tmpj < n; tmpj++)
 //          {
 //                printf("%f  ", V[tmpi*n + tmpj]);
 //          }
 //          printf("\n");
 //    }
 //    printf("\n\n");
       
                   // Implicit QL transformation.
       
                   p = d[m];
                   //Real c = 1.0;
                   //Real c2 = c;
                   //Real c3 = c;
                   //Real el1 = e[l+1];
                   //Real s = 0.0;
                   //Real s2 = 0.0;
                   c = 1.0;
                   c2 = c;
                   c3 = c;
                   el1 = e[l+1];
                   s = 0.0;
                   s2 = 0.0;
                   
                   
                   for (i = m-1; i >= l; i--)
                   {
   //                    printf ("p top: %E\n", p);
                       c3 = c2;
                       c2 = c;
                       s2 = s;
                       g = c * e[i];
                       h = c * p;
                       r = hypot(p,e[i]);
     //                  printf("r = %E\n", r);
                       e[i+1] = s * r;
                       s = e[i] / r;
                       c = p / r;
                       p = c * d[i] - s * g;
       //                printf ("d[i]: %E\n", d[i]);
        //               printf ("p bot: %E\n", p);
         //              printf ("c*d[i]: %E\n", c * d[i]);
        //               printf ("s*g: %E\n",s * g);
        //               printf ("s: %E\n", s);
        //               printf ("c: %E\n", c);
        //               printf ("g: %E\n", g);
                       d[i+1] = h + s * (c * g + s * d[i]);
                       
                       
                       
    //                   printf("d3.5: %E %E %E %E\n",d[0],d[1],d[2],d[3]);
   //  printf("e3.5: %E %E %E %E\n",e[0],e[1],e[2],e[3]);
  //   printf("V3.5:\n");
//          for(tmpi = 0; tmpi < n; tmpi++)
//     {
//           for(tmpj = 0; tmpj < n; tmpj++)
//           {
//                 printf("%f  ", V[tmpi*n + tmpj]);
//           }
//           printf("\n");
//     }
//     printf("\n\n");
                       
                       
                       
                       
                       // Accumulate transformation.
                       
                       for (k = 0; k < n; k++)
                       {
                           //h = V[k][i+1];
                           h = V[k*n +i+1];
                           
                           //V[k][i+1] = s * V[k][i] + c * h;
                           V[k*n +i+1] = s * V[k*n +i] + c * h;
                           
                           //V[k][i] = c * V[k][i] - s * h;
                           V[k*n +i] = c * V[k*n +i] - s * h;
                       }
                       
                       
 //                      printf("d3.75: %E %E %E %E\n",d[0],d[1],d[2],d[3]);
 //    printf("e3.75: %E %E %E %E\n",e[0],e[1],e[2],e[3]);
  //   printf("V3.75:\n");
  //        for(tmpi = 0; tmpi < n; tmpi++)
  //   {
  //         for(tmpj = 0; tmpj < n; tmpj++)
  //         {
  //               printf("%f  ", V[tmpi*n + tmpj]);
  //         }
  //         printf("\n");
  //   }
  //   printf("\n\n");
                   }
                   
                   p = -s * s2 * c3 * el1 * e[l] / dl1;
                   e[l] = s * p;
                   d[l] = c * p;
                   
                   // Check for convergence.
   
            }
            while (fabs(e[l]) > eps*tst1);
         }
         d[l] = d[l] + f;
         e[l] = 0.0;
   }
   
 //  printf("d4: %E %E %E %E\n",d[0],d[1],d[2],d[3]);
 //    printf("e4: %E %E %E %E\n",e[0],e[1],e[2],e[3]);
 //    printf("V4:\n");
 //         for(tmpi = 0; tmpi < n; tmpi++)
  //   {
  //         for(tmpj = 0; tmpj < n; tmpj++)
  //         {
  //               printf("%f  ", V[tmpi*n + tmpj]);
   //        }
   //        printf("\n");
  //   }
  //   printf("\n\n");
     
      // Sort eigenvalues and corresponding vectors.
//   This should be unnecessary for us, but it'll make testing easier
      for (i = 0; i < n-1; i++) {
         int k = i;
         Real p = d[i];
         for (j = i+1; j < n; j++) {
            if (d[j] < p) {
               k = j;
               p = d[j];
            }
         }
         if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 0; j < n; j++) {
               //p = V[j][i];
               p = V[j*n + i];
               
               //V[j][i] = V[j][k];
               V[j*n + i] = V[j*n + k];
               
               V[j*n + k] = p;
               V[j*n + k] = p;
            }
         }
      }
      
}




//----------------------My Original C Cross Section Code------------------------

/*
I may want to change this in the future, but to implement 2D arrays here, I just
used arrays of pointers, so that I could index like M[i][j].
*/


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

int main()
{
    
}
