/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation. 
 * Any use, reproduction, disclosure, or distribution of this software 
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA) 
 * associated with this source code for terms and conditions that govern 
 * your use of this NVIDIA software.
 * 
 */

#include "../common/book.h"

__global__ void tred2(Real *V, Real *d, Real *e, int n);
__global__ void tql2(Real *V, Real *d, Real *e, int n);

typdef float Real;

//-----------------------Eigenvalue/vector Code---------------------------------

/*I am taking the releveant code from TNT and adapting it to our problem in C.
I am taking only the code necessary to find the eigenvectors of a hermitian
matrix.  The TNT code is only for real matrices, but luckilly a Hermitian
matrix can be converted to a real matrix to find the eigenvalues/vectors.
*/

void test_func()
{
     int i,j;
     int n = 4;//The TNT code uses this variable name
     Real *dev_V, *dev_d, *dev_e;
     
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
     
     
     //hardcode my real matrix here for testing
     n = 4;
     V = (Real *) malloc(sizeof(Real)*n*n);
     V[0] = 1;
     V[1] = 5;
     V[2] = 3;
     V[3] = 7;
     V[4] = 5;
     V[5] = 4;
     V[6] = 8;
     V[7] = 3;
     V[8] = 3;
     V[9] = 8;
     V[10] = 1;
     V[11] = 7;
     V[12] = 7;
     V[13] = 3;
     V[14] = 7;
     V[15] = 5;
     
     
     printf("The real Matrix:\n");
     for(i = 0; i < n; i++)
     {
           for(j = 0; j < n; j++)
           {
                 printf("%f  ", V[i*n + j]);
           }
           printf("\n");
     }
     
     //Allocate Space on the GPU
     HANDLE_ERROR( cudaMalloc( (void**)&dev_V, n * n * sizeof(Real) ) );
     HANDLE_ERROR( cudaMalloc( (void**)&dev_d, n * sizeof(Real) ) );
     HANDLE_ERROR( cudaMalloc( (void**)&dev_e, n * sizeof(Real) ) );
     
     //Copy V over to the GPU
     HANDLE_ERROR( cudaMemcpy( dev_V, V, n*n * sizeof(Real), cudaMemcpyHostToDevice ) );
     
     //Tridiagonalize (on GPU)
     tred2<<<1,1>>>(dev_V, dev_d, dev_e, n);
     
     //Diagonalize (on GPU)
     tql2<<<1,1>>>(dev_V, dev_d, dev_e, n);
     
     //Copy V,d, and e
    HANDLE_ERROR( cudaMemcpy( V, dev_V, n*n * sizeof(Real),
                              cudaMemcpyDeviceToHost ) );
    HANDLE_ERROR( cudaMemcpy( d, dev_d, n*n * sizeof(Real),
                              cudaMemcpyDeviceToHost ) );
    HANDLE_ERROR( cudaMemcpy( e, dev_e, n*n * sizeof(Real),
                              cudaMemcpyDeviceToHost ) );
     
     //test
     printf("\n\neigs:\n");
     for(i = 0; i < n; i++)
     {
           //printf("eval: %f + %fi\n", d[i], e[i]);
           //printf("evec row:\n");
           for(j = 0; j < n; j++)
           {
                 printf("%f ", V[i*n + j]);
           }
           printf("\n");
     }
     for(i = 0; i < n; i++)
     {
           printf("eval: %f + %fi\n", d[i], e[i]);
     }
     
     
     free(V);
     free(d);
     free(e);
     
    // free the memory allocated on the GPU
    HANDLE_ERROR( cudaFree( dev_V ) );
    HANDLE_ERROR( cudaFree( dev_d ) );
    HANDLE_ERROR( cudaFree( dev_e ) );
}



//Adapted from the TNT code (Eigenvalue class) in jama_eig.h
__global__ void tred2(Real *V, Real *d, Real *e, int n)
{
//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

    int i,j,k;
    Real scale, h, f, g, hh;
    
    for (j = 0; j < n; j++)
    {
        //d[j] = V[n-1][j];
        d[j] = V[(n-1)*n+j];
    }

    // Householder reduction to tridiagonal form.

    for (i = n-1; i > 0; i--)
    {

        // Scale to avoid under/overflow.
        
        //Real scale = 0.0;
        //Real h = 0.0;
        
        //I haven't looked through all the code for all the uses of these, but
        //I think it's very unlikely that they can't be overwritten for each
        //iteration.
        scale = 0.0;
        h = 0.0;
        
        for (k = 0; k < i; k++)
        {
            scale = scale + fabs(d[k]);
        }
        if (scale == 0.0)
        {
           
           
            e[i] = d[i-1];
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
        
            // Generate Householder vector.
            for (k = 0; k < i; k++)
            {
                d[k] /= scale;
                h += d[k] * d[k];
            }

            //Real f = d[i-1];
            //Real g = sqrt(h);
            f = d[i-1];
            g = sqrt(h);


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
            
            
            
            f = 0.0;
            
            for (j = 0; j < i; j++)
            {
                e[j] /= h;
                f += e[j] * d[j];
            }
            //Real hh = f / (h + h);
            hh = f / (h + h);
//            printf("\nhh: %E\n", hh);
            
            for (j = 0; j < i; j++)
            {
                e[j] -= hh * d[j];
            }
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
        }
        d[i] = h;
    
    
    }
    
    
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
    
    
}




//Adapted from function in Eigenvalue class in TNT code

// Symmetric tridiagonal QL algorithm.
   
__global__ void tql2(Real *V, Real *d, Real *e, int n)
{
     //  This is derived from the Algol procedures tql2, by
     //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
     //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
     //  Fortran subroutine in EISPACK.
     

     int i, l, m, iter, k, j;
     Real f, tst1, eps, g, p, r, dl1, h, c, c2, c3, el1, s, s2;
     
     
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
                   
       
                   // Implicit QL transformation.
       
                   p = d[m];
                   c = 1.0;
                   c2 = c;
                   c3 = c;
                   el1 = e[l+1];
                   s = 0.0;
                   s2 = 0.0;
                   
                   
                   for (i = m-1; i >= l; i--)
                   {
                       c3 = c2;
                       c2 = c;
                       s2 = s;
                       g = c * e[i];
                       h = c * p;
                       r = hypot(p,e[i]);
                       e[i+1] = s * r;
                       s = e[i] / r;
                       c = p / r;
                       p = c * d[i] - s * g;
                       d[i+1] = h + s * (c * g + s * d[i]);
                       
                       
                       
                       
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



int main( void ) {
    test_func();

    return 0;
}
