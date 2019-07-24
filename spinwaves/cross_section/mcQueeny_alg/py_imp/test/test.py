import numpy as np
import cmath


sigmas = [1,-1]
q = 2*np.pi/3.0

M = np.array([
[2, -1-np.exp(-1j*q)],
[1+np.exp(1j*q), -2]])
#M = np.array([
#[1,2,3 + 1j],
#[2 + 1j, -1+3j, 1-1j],
#[3,4+1j, 1-1j]])


w, evec = np.linalg.eig(M)

print "w: ", w
print "evec: ", evec



#Get number of atoms
numAtoms = 2

T = np.empty((numAtoms, numAtoms), dtype = complex)#normalized eigenvectors

#Normalize evec to get T_ni
for n in range(0, numAtoms):
    sgn_ln = w[n]/abs(w[n])
    tmp = 0.0
    for i in range(0, numAtoms):
        tmp += sigmas[i] * abs(evec[i][n])**2
    C = sgn_ln/tmp
    #print "C = ", sgn_ln, "/", tmp, " = ", C
    for i in range(0, numAtoms):
        T[i][n] = evec[i][n]*cmath.sqrt(C)
        #print "type(evec[i][n]): ", type(evec[i][n])
        #print "T[",i,"][",n,"] = ", evec[i][n], " * ", cmath.sqrt(C), " = " , T[i][n]
        
print "T: ", T

print "T1: ", sigmas[0]*(abs(T[0][0])**2) + sigmas[1]*(abs(T[1][0])**2)# + sigmas[2]*abs(T[2][1])**2

#Check which way the indices go on the eigenvectors
#V1 = []
#V1.append(T[0][0])
#V1.append(T[1][0])
#V1 = np.array(V1).transpose()
#V2 = []
#V2.append(T[0][1])
#V2.append(T[1][1])
#V2 = np.array(V2).transpose()
V1 = np.array([
[T[0][0]],
[T[1][0]]])
V2 = np.array([
[T[0][1]],
[T[1][1]]])
print "V1: ", V1
print "V2: ", V2
print "np.dot(M,V1) - w[0]*V1 : ", np.dot(M,V1) - w[0]*V1
print "np.dot(M,V2) - w[1]*V2 : ", np.dot(M,V2) - w[1]*V2


v1 = np.array([evec[0][0], evec[1][0]])
print np.dot(M,v1.transpose())
print w[0]*v1




Q = np.array([0,1.0,1.0])
spin = np.array([0,0,5])
mu_hat = spin/np.sqrt(np.dot(spin,spin))
polarization_factor = 0.5*(1.0 + (np.dot(mu_hat,Q)**2)/np.dot(Q,Q))
print "mu_hat: ", mu_hat
print polarization_factor

