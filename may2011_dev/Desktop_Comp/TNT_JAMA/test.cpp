#include "jama_eig.h"

//My first C++ program (since 10th grade which didn't count)

//Testing the eigenvector calculation in jama_eig.h

//my main function
int main()
{
    TNT::Array2D<double> testAr;
    //testAr = Array2D<double>(3,3);
    testAr = Array2D<double>(4,4);
//    testAr[0][0] = 1.0;
//    testAr[0][1] = 1.0;
//    testAr[0][2] = 0.0;
//    testAr[0][3] = -1.1;
//    
//    testAr[1][0] = 1.0;
//    testAr[1][1] = 0.3;
//    testAr[1][2] = 1.1;
//    testAr[1][3] = 0.0;
//    
//    testAr[2][0] = 0.0;
//    testAr[2][1] = 1.1;
//    testAr[2][2] = 1.0;
//    testAr[2][3] = 1.0;
//    
//    testAr[3][0] = -1.1;
//    testAr[3][1] = 0.0;
//    testAr[3][2] = 1.0;
//    testAr[3][3] = 0.3;

testAr[0][0] = 1.0;
    testAr[0][1] = 2.1;
    testAr[0][2] = 0.0;
    testAr[0][3] = -0.4;
    
    testAr[1][0] = 2.1;
    testAr[1][1] = 3.2;
    testAr[1][2] = 0.4;
    testAr[1][3] = 0.0;
    
    testAr[2][0] = 0.0;
    testAr[2][1] = 0.4;
    testAr[2][2] = 1.0;
    testAr[2][3] = 2.1;
    
    testAr[3][0] = -0.4;
    testAr[3][1] = 0.0;
    testAr[3][2] = 2.1;
    testAr[3][3] = 3.2;
    
    //JAMA::Eigenvalue<double> eig;
    JAMA::Eigenvalue<double>eig = JAMA::Eigenvalue<double>::Eigenvalue(testAr);

    //eig = Eigenvalue<double>(testAr);;
    //eig = Eigenvalue(testAr);
    cout << testAr[0][0] << " " << testAr[0][1] << " " << testAr[0][2] << " " << testAr[0][3] << "\n";
    cout << testAr[1][0] << " " << testAr[1][1] << " " << testAr[1][2] << " " << testAr[1][3] << "\n";
    cout << testAr[2][0] << " " << testAr[2][1] << " " << testAr[2][2] << " " << testAr[2][3] << "\n";
    cout << testAr[3][0] << " " << testAr[3][1] << " " << testAr[3][2] << " " << testAr[3][3] << "\n";
    cout << "running\n";
    eig.getV(testAr);
    cout << testAr[0][0] << " " << testAr[0][1] << " " << testAr[0][2] << " " << testAr[0][3] << "\n";
    cout << testAr[1][0] << " " << testAr[1][1] << " " << testAr[1][2] << " " << testAr[1][3] << "\n";
    cout << testAr[2][0] << " " << testAr[2][1] << " " << testAr[2][2] << " " << testAr[2][3] << "\n";
    cout << testAr[3][0] << " " << testAr[3][1] << " " << testAr[3][2] << " " << testAr[3][3] << "\n";
    cout << "eigvals:\n";
    
    TNT::Array1D<double> d = Array1D<double>(3);
    TNT::Array1D<double> e = Array1D<double>(3);
    
    eig.getRealEigenvalues(d);
    eig.getImagEigenvalues(e);
    
    cout << d[0] << " + " << e[0] << "j\n"; 
    cout << d[1] << " + " << e[1] << "j\n";
    cout << d[2] << " + " << e[2] << "j\n";
    cout << d[3] << " + " << e[3] << "j\n";

    
    //cout << "Hello world!";
    //So it looks like what I need to do is instantiate the Eingenvalue class
    //with a 2D array instance and it will automatically calculate the eigenvector
    //matrix in the constructor, which I can then get with getV()
    
    //So first, make a 2D array
    //Unfortunately, it looks like this only supports real matrices.
    
    return 0;
}
