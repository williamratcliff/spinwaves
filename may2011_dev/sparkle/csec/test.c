#include <stdio.h>

void test()
{
     printf("success!\n");
}


//Adapted from the TNT code (Eigenvalue class) in jama_eig.h
__global__ void tred2(Real *V, Real *d, Real *e, int n)
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
