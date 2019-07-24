import numpy as np
cimport numpy as np
import time
from time import clock
from spinwaves.spinwavecalc.readfiles import atom, readFiles
from periodictable import elements

cimport cMcQueenyAlg as calg

cdef int* intListToCArray(list):
    cdef int* array = calg.allocateIntArray(len(list))
    for i in range(len(list)):
        array[i] = list[i]
    return array

cdef double* dblListToCArray(list):
    cdef double* array = calg.allocateDoubleArray(len(list))
    for i in range(len(list)):
        array[i] = list[i]
    return array

#takes a list of python SimpleAtoms and returns a list in C
cdef calg.Atom* createCAtomList(list, jnums, jmats, spinAxis):
    cdef calg.Atom *clist = calg.createAtomList(len(list))
    #cdef double *dblJs
    for i in range(len(list)):
        dblJs = []
        for j in range(len(list[i].interactions)):
            dblJs.append(get_Jij(jnums, jmats, list[i].interactions[j])[1][1])
        calg.addToList(clist, i, len(list[i].neighbors), intListToCArray(list[i].neighbors), dblListToCArray(dblJs), list[i].spin[spinAxis], list[i].pos[0], list[i].pos[1], list[i].pos[2])
        #dblJs = calg.allocateDoubleArray(len(list[i].interactions))
        #for j in range(len(list[i].interactions)):
        #    print "here"
        #    dblJs[j] = get_Jij(jnums, jmats, list[i].interactions[j])
        #    print "here1"
        #calg.addToList(clist, i, len(list[i].neighbors), intListToCArray(list[i].neighbors), dblJs, list[i].spin[spinAxis], list[i].pos[0], list[i].pos[1], list[i].pos[2])
    return clist


cdef np.ndarray[np.complex128_t, cast=True] sqr_mat_to_np(calg.SquareMatrix *sqMat):
    cdef np.ndarray[np.complex128_t, ndim=2] a = np.zeros(((sqMat[0]).dim, (sqMat[0]).dim), dtype=np.complex128)
    cdef int i,j
    for i in range((sqMat[0]).dim):
        for j in range((sqMat[0]).dim):
            a[i][j] = (sqMat[0]).mat[i][j].real + 1j*(sqMat[0]).mat[i][j].imaginary
    return a

def get_Jij(jnums, jmats, num):
    """Returns the 3x3 numpy array corresponding to the matrix number, num from an atoms list
    of interactions."""
    index = jnums.index(num)
    return np.float64(np.array(jmats[index])) #This would be a sympy Matrix if I didn't use the np wrapper

cdef calg.ComplexNum** convert_evec(evec):
    """Converts the numpy array of eigenvecttors to a 2D C array of calg.ComplexNums (ComplexNum **)."""
    cdef calg.ComplexNum **cEvec = calg.allocate_evec(len(evec))
    for i in range(len(evec)):
        for j in range(len(evec)):
            cEvec[i][j].real = evec[i][j].real
            cEvec[i][j].imaginary = evec[i][j].imag
    return cEvec

cdef double* convert_w(w):
    for val in w:
        if val.imag > (val.real*1e-10):
            raise Exception("Complex Eigenvalue Error!")
    return dblListToCArray(w)


#Old code
def mcQueeny_case(interactionfile, spinfile, tau_list, temperature,
                      direction={'kx':1,'ky':0,'kz':0}, hkl_interval=[1e-3,2*np.pi,100],
                      omega_interval=[0,5,100], eig_eps = 0.001):

    print "starting McQueeny Case"
    initial_time = time.clock()

    #to plot eigenvector components
#    q1 = []
#    q2 = []

    #to plot w's
    #w1_list = []
    #w2_list = []
    #w3_list = []
    #w4_list = []
    #w5_list = []
    #w6_list = []
    #w7_list = []
    #w8_list = []
#    w_output_list = []

    #plot the structure factor to compare with lovesey (In McQueeny Paper)
#    strfac_sums = []
#    strfac1 = []
#    strfac2 = []

    atom_list, jnums, jmats, N_atoms_uc, numCutOff, magCellSize = readFiles(interactionfile,spinfile, rtn_cutOffInfo = True)

    print "creating the list in C"
    #C implementation
    cdef calg.Atom *clist = createCAtomList(atom_list, jnums, jmats, 2)

    h_list, k_list, l_list = generate_hkl(hkl_interval, direction)
    e = generate_wt(omega_interval)
    ff_list = generate_form_factors(N_atoms_uc, atom_list, hkl_interval)
    s = []
    basis = []
    for a in atom_list:
        s.append(a.spin[2])#z component of spin (assume aligned with z)
        #basis.append(a.pos)
        basis.append(a.pos - magCellSize) #translate  back to 0,0,0
        #normalize? NO!
        #basis.append([a.pos[0] - int(a.pos[0]), a.pos[1]-int(a.pos[1]), a.pos[2] - int(a.pos[2])])
    #print "\n\nbasis vectors:\n"

    #mess with basis vectors and spin
    #basis[0] = [0.25,0,0]
    #basis[0] = [1.0,0,0]
    #tmp = s[0]
    #s[0] = s[1]
    #s[1] = tmp
    #for pos in basis:
    #    print pos
    #print "\n\nspins:\n"
    #for spin in s:
    #    print spin
    mesh = np.zeros((len(h_list), len(e)))

    #output q, w, and eigenvectors before and after normalization to compare with McQueeny's
#    f = open("/home/tom/Desktop/Python_Eigs.txt", 'w')
#    f.write("q w eigVec_before_norm T\n\n")

#    f2 = open("/home/tom/Desktop/Values.txt", 'w')
#    f2.write("Q w wvec (wvec norm) e^(-i*Q.d_i)\n\n")


    #create the T_ni matricies (eigenvector matrices)
    #This algorithm uses a differnt matrix:
    #evecs, mats = calc_eig_mats_numerically(Hsave, h_list,k_list,l_list)
    #eigvals = calc_eigs_numerically(Hsave, h_list, k_list, l_list)#This is innefficient, eigs already calculated and thrown away in calc_eig_mats_..
    #w_list = gen_eigs_with_sign(Hsave, h_list, k_list, l_list, 0.001)

    #declare all the C variables I will need (that aren't primitives)
    cdef calg.SquareMatrix *m
    cdef np.ndarray[np.complex128_t, ndim=2] M
    cdef double **cMesh = calg.allocate_mesh(len(h_list), len(e))
    #cdef calg.ComplexNum **cEvec
    #cdef double *cW
    #cdef double *cE = dblListToCArray(e)

    for Qindex in range(len(h_list)):
        #print "calculating Q point..."
        Q = (np.array([h_list,k_list,l_list]).transpose())[Qindex]
#        f2.write(str(Q))
        #CALL SPINWAVE(q,w,evec) - get w, evec for given q
        #M = create_matrix(atom_list, jnums, jmats, numCutOff, magCellSize, Q-np.array(tau_list))

        #C implementation
        #print "creating matrix..."
        m = calg.createMatrix(clist,numCutOff, Q[0],Q[1],Q[2], magCellSize[0], magCellSize[1], magCellSize[2])
        #print "converting matrix to numpy matrix"
        M = sqr_mat_to_np(m)

        #print "finding eigs..."
        w, evec = np.linalg.eig(M)



#        f2.write(" " + str(w))
#        f2.write(" " + str(evec.flatten()))

        #                evec = evecs[Qindex]
        #                w = w_list[Qindex]




      #  print "\nw = ", w
        #w1_list.append(w[0])
        #w2_list.append(w[1])
        #w3_list.append(w[1])
        #w4_list.append(w[1])
        #w5_list.append(w[1])
        #w6_list.append(w[1])
        #w7_list.append(w[1])
        #w8_list.append(w[1])
#        w_output_list.append(w)


        ff = ff_list[Qindex]

        calg.calc_cs(cMesh, Qindex, convert_evec(evec), convert_w(w), numCutOff, clist, dblListToCArray(e), len(e), ff, Q[0], Q[1], Q[2], magCellSize[0], magCellSize[1], magCellSize[2], LIFETIME_VALUE)

##        #Get number of atoms
##        numAtoms = numCutOff
##        #numAtoms = np.sqrt(evec.size)#This should be the number of atoms that I need to loop through.
##        #if int(numAtoms) != numAtoms:#I beleive t_mat should always be square, but lets check
##        #    raise Exception("t_mat not square!")
##
##        #Normalize evec to get T_ni
##        #for i in range(0, numAtoms):
##        #    normi = 0
##        #    for j in range(0, numAtoms):
##        #        sigma_j = s[j]/np.abs(s[j])
##        #        normi += sigma_j * np.abs(evec[j][i])**2
##        #    for j in range(0, numAtoms):
##        #        evec[i][j] = evec[j][i]/np.sqrt(np.abs(normi))
##
##        evec2 = evec.copy()
##        #I'm going to write this from scratch
##      #  print "evec before = ", evec
##        T = evec
##        sigma = []
##        for i in range(numAtoms):
##            sigma.append(s[i]/np.abs(s[i]))
##        for n in range(numAtoms):
##            norm = 0
##            for i in range(numAtoms):
##                norm += sigma[i]*(np.abs(evec2[n][i])**2)
##            for i in range(numAtoms):
##      #          print "norm : ", norm
##      #          print "w sign: ", np.sign(w[n])
##                #T[n][i] = T[n][i]*np.sqrt((w[n]/np.abs(w[n]))/norm + 0j)
##                T[n][i] = evec2[n][i]/np.sqrt(np.abs(norm))
##        evec = T
###        f2.write(" " + str(T.flatten()))
##        #evec = evec2
###        q1.append(evec[0][0])
###        q2.append(evec[0][0])
##
##        #output eigs to file
###        f.write(str(Q))
###        f.write(" w:" + str(w))
###        f.write("\nmat: " + str(M))#the matrix for which the eigs were found
###        f.write("\nbefore_norm: " + str(evec2.flatten()))
###        f.write("\nafterNorm: " + str(evec.flatten()) + "\n\n")
##
##        #evec = evec.transpose()
##        #test = 0
##        #for i in range(numAtoms):
##        #    sigma_i = s[j]/np.abs(s[j])
##        #    test += sigma_i * np.dot(evec[i],evec[i])
##        #    print "sigma_i = ", sigma_i
##        #print "\n\nsum{sigma_i * |T_ni|^2} = ", test, "\n\n"
##
##        #                print "evec = " , evec
##        #                print "sigma1: ", s[0]/np.abs(s[0])
##        #                print "sigma2: ", s[1]/np.abs(s[1])
##        #                print "sign lambda_0: ", (s[0]/np.abs(s[0]))*(evec[0][0]**2) + (s[1]/np.abs(s[1]))*(evec[0][1]**2)
##        #                print "sign lambda_1: ", (s[0]/np.abs(s[0]))*(evec[1][0]**2) + (s[1]/np.abs(s[1]))*(evec[1][1]**2)
##
##
##        #using only positive eigenvalue
##        #w = [w[0]]
##        #evec = [evec[0]]
##        #print "eigenvalues used: ", w
##        #print "eigenvectors used: ", evec
##
##
##
##        strfac = np.zeros(numAtoms)
##
##        #for testing:
##        ff = 1
##
##        for i in range(numAtoms):
##            dum = 0
##            for j in range(numAtoms):
##                dum += (s[j]/np.sqrt(np.abs(s[j]))) * ff * evec[j][i] #* np.exp(1.0j*
##                        #(-Q[0]*basis[j][0] -Q[1]*basis[j][1] -Q[2]*basis[j][2]))
###                print "\n\n\nwritting to f2:\n\n\n", np.exp(1.0j*
###                        (-Q[0]*basis[j][0] -Q[1]*basis[j][1] -Q[2]*basis[j][2]))
###                print "j: ", j
###                print "basis: ", basis[j]
###                print "Q.di", (-Q[0]*basis[j][0] -Q[1]*basis[j][1] -Q[2]*basis[j][2])
##                exp = np.exp(1.0j*(-Q[0]*basis[j][0] -Q[1]*basis[j][1] -Q[2]*basis[j][2]))
###                f2.write(" " + str(exp))
##            strfac[i] = np.abs(dum)**2
##       #     print "Q: ", Q, "    strfac[",i,"]: ", strfac[i]
###            f2.write("\n\n")
###        strfac_sums.append(0.0);
###        strfac1.append(strfac[0])
###        strfac2.append(strfac[1])
##        sqw = np.zeros(len(e))
##        for j in range(numAtoms):
###            strfac_sums[Qindex] +=strfac[j]
##            for i in range(len(e)):
##                y = LIFETIME_VALUE/2
##                lorentzian = (1/np.pi)*y/((e[i]- w[j])**2+y**2)
##                #test
##                #strfac[j] = ff*1
##                sqw[i] = sqw[i] + strfac[j] * lorentzian #he has some 'bose' factor in here
##                #more test
##                #sqw[i] = sqw[i] + strfac[j]
##        mesh[Qindex] = sqw #* 0.5*(1 + Q[0]**2) #along x, Q_hat is just 1?

    #f.close()
    #f = open("/home/tom/Desktop/output.txt", 'w')
    #for q in range(0,len(h_list)):
        #for w in range(0,len(e)):
            #f.write(str(h_list[q]))
            #f.write(" ")
            #f.write(str(e[w]))
            #f.write(" ")
            #f.write(str(mesh[q][w]))
            #f.write("\n")
    #f.close()
    #f = open("/home/tom/Desktop/output2.txt", 'w')
    #for q in range(0,len(h_list)):
            #f.write(str(h_list[q]))
            #f.write(" ")
            #f.write(str(np.abs(q1[q])))
            #f.write(" ")
            #f.write(str(np.abs(q2[q])))
            #f.write("\n")
    #f.close()
    #f = open("/home/tom/Desktop/output3.txt", 'w')
    #for q in range(0,len(h_list)):
            #f.write(str(h_list[q]))
            #f.write(" ")
            #f.write(str(np.real(w1_list[q])))
            #f.write(" ")
            #f.write(str(np.real(w2_list[q])))
            #f.write(" ")
            #f.write(str(np.real(w3_list[q])))
            #f.write(" ")
            #f.write(str(np.real(w4_list[q])))
            #f.write(" ")
            #f.write(str(np.real(w5_list[q])))
            #f.write(" ")
            #f.write(str(np.real(w6_list[q])))
            #f.write(" ")
            #f.write(str(np.real(w7_list[q])))
            #f.write(" ")
            #f.write(str(np.real(w8_list[q])))
            #f.write("\n")
    #f.close()
    #f = open("/home/tom/Desktop/strfacs.txt", 'w')
    #for q in range(0,len(h_list)):
            #f.write(str(h_list[q]))
            #f.write(" ")
            #f.write(str(k_list[q]))
            #f.write(" ")
            #f.write(str(l_list[q]))
            #f.write(" ")
            #f.write(str(strfac1[q]))
            #f.write(" ")
            #f.write(str(strfac2[q]))
            #f.write(" ")
            #f.write(str(strfac_sums[q]))
            #f.write("\n")
    #f.close()
    #plt.plot(w1_list)
    #plt.plot(w2_list)
    #plt.show()
    print "Run Time: ", time.clock()-initial_time, " seconds"



#The following functions were copied from csection_calc, but also left in that file becuase they are useed by both
#the MyQueeny algorithm and ours.  I duplicated the code becuase the code in this file is intended to be converted to Cython
#and I assume that talking between cython and pthon is slow.




#A copy of this was also left in cesction-calc because I assume talking between python and cyton is slow
def generate_hkl(hkl_interval, direction):

    #TODO - THIS DOES NOT WORK FOR ALL CASES! only when directional components are the same or 0
    #Fixed - testing

    #What if other points are 0 (not the first or last points)?
    #For that to happen they would need to generate an odd number of points across an interval -x to x
    if hkl_interval[0] == 0:
        hkl_interval[0] = 1e-3
    elif hkl_interval[1] == 0:
        hkl_interval[1] = -1e-3

    print "hkl_interval: ", hkl_interval
    print "direction: ", direction

    #We should be multipliying these linspaces by the components, kx,ky,kz, where those are
    #normalized so that kx**2 + ky**2 + kz**2 = 1 (assuming that h**2+k**2+l**2=q**2), I'm not
    #entirely sure if this works in inverse space.
    norm = np.sqrt(direction['kx']**2 + direction['ky']**2 + direction['kz']**2)

    if direction['kx']:
        h_list = np.linspace(hkl_interval[0],hkl_interval[1],hkl_interval[2])*direction['kx']/norm
    else:
        h_list = np.zeros(hkl_interval[2])
    if direction['ky']:
        k_list = np.linspace(hkl_interval[0],hkl_interval[1],hkl_interval[2])*direction['ky']/norm
    else:
        k_list = np.zeros(hkl_interval[2])
    if direction['kz']:
        l_list = np.linspace(hkl_interval[0],hkl_interval[1],hkl_interval[2])*direction['kz']/norm
    else:
        l_list = np.zeros(hkl_interval[2])

    #print "h,k,l points: "
    #for i in range(len(h_list)):
    #    print "(",h_list[i],", ", k_list[i], ", ", l_list[i], ")"

    return h_list, k_list, l_list

#A copy of this was also left in cesction-calc because I assume talking between python and cyton is slow
def generate_wt(omega_interval):

    return np.linspace(omega_interval[0],omega_interval[1],omega_interval[2])

#A copy of this was also left in cesction-calc because I assume talking between python and cyton is slow
def generate_form_factors(N_atoms_uc, atom_list, hkl_interval):
    # Form Factor
    eval_pnts = np.linspace(hkl_interval[0],hkl_interval[1],hkl_interval[2])
    ff_list = []
    for i in range(N_atoms_uc):
        try:
            elSym = atom_list[i].label
            elSym = elements.symbol(elSym)
            elMass = atom_list[i].massNum
            if elMass == None: #The mass number will now generally be None -Tom 1/28/11
                el = elSym
            else:
                el = elSym[elMass]
            print el
            val = atom_list[i].valence
            if val != None:
                Mq = el.magnetic_ff[val].M_Q(eval_pnts)
            else:
                Mq = el.magnetic_ff[0].M_Q(eval_pnts)
        except:
            Mq = np.zeros(len(eval_pnts))
        #Test to see how things look without ff
        #Mq = np.ones(len(eval_pnts))

        ff_list = Mq
#    fake = np.ones(len(kaprange))
    return np.array(ff_list)

from csection_calc import LIFETIME_VALUE