"""This is an adaptation of McQueeny's algorithm to python.  Right now, the
spins must be aligned along the z axis.  I plan on allowng them to be aligned
along any axis by rotating my whole coordinate system so that the z axis is
defined to be along the direction of the spins.  I am also setting the form
factors to one for now.

This technique only works for lattices with collinear spins.  Therefore, the
interaction matrices are assumed to be diagonal, representing either
ferromagnetic or anti-ferromagnetic interactions.

For now, I am assuming the Bose occupation factor, n_qn in the paper, is zero.
In McQueeny's code he has:  bose = 1. + 1./(EXP(11.6*w(j)/temp) - 1.)

This is (1 + )the Bose-Einstein distribution for identical bosons
(pg. 241 Griffiths) where 11.6 = 1/kb where kb is in units of meV/K.

Once I make this support spin directions other than the z axis, I should
probably add a function to test if all the spins are aligned within some
tolerance and then pick the axis direction - I will need to change the method
with which I find the sigmas.  I should be able to either find a rotation
matrix and rotate all vector inputs to teh algorithm before running it so that
the spins align with the z axis, or leave all the vectors the same, but change
the sigma calculation algorithm and change the polarization factors.

1/13/12  The structure factor results of this ( S(q), excluding delta functions,
r_0, k'/k, formfactors...) match what I get for the simple anti-ferromagnetic
chain - in file simple_anti-ferro_test.py.

Tom
1/4/12"""

import cmath
import time
import numpy as np
from matplotlib import pyplot as plt

import sys
sys.path.append(r"/home/tom/Desktop/spinwaves/trunk")
from spinwaves.spinwavecalc.readfiles import atom, readFiles
from periodictable import elements

#just a test
def is_real(M):
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.abs(M[i][j].real) < np.abs(M[i][j].imag)*1000:
                print M[i][j]
                return False
    return True

def is_herm(M):
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            assert(M[i][j].real == M[j][i].real and M[i][j].imag == -M[j][i].imag)
    print "hermitian!"
    #time.sleep(1)


def dbl_equal(d1, d2, maxAbsoluteError = 1.0e-12, maxRelativeError = 1.0e-7):
    """A good discussion of floating point equality comparison can be found here:
    http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
    
    I adapted this form one of his C AlmostEqual functions."""
    if abs(d1-d2) < maxAbsoluteError:
        return True
    if abs(d2) > abs(d1):
        relativeError = abs((d1 - d2) / d2)
    else:
        relativeError = abs((d1 - d2) / d1)
    if relativeError <= maxRelativeError:
        return True; 
    return False
    


def csection(interactionfile, spinfile, tau_list, temp = 0.0001, HWHM = 1/np.pi,
            direction={'kx':1,'ky':0,'kz':0}, hkl_interval=[1e-3,2*np.pi,100],
            omega_interval=[0,5,100], kb = 0.086173423):
    """
    kb is the boltzman constant in whatever units the user is using.  The
    default is 0.086173423 meV/K.
    
    HWHM is the half width at half maximum of the lorentzian broadening
    function (the resolution of the energy)."""

    initial_time = time.clock()
    

    atom_list, jnums, jmats, N_atoms_uc, numCutOff, magCellSize = readFiles(interactionfile,spinfile, rtn_cutOffInfo = True)
    #Atoms in the first interaction cell come first in the list (atom_list),so
    #I can iterate through the first numCutOff atoms in the list top iterate
    #only through the atoms in the first interaction cell - I beleive 1/4/12
    
    #Right now (1/12/12) a major issue is that this method requires the size
    #of the magnetic cell and the cutOffCell size is not necesssarilly the
    #same as the magenetic cell size.  In fact, even in the very simple
    #ferromagnetic chain it is not. - Same issue with numCutOff - should be
    #numMagCell
    
    #Simple Ferro Case
    #magCellSize = np.array([1,1,1])
    #numCutOff =1
    
    
    #testing
    #numCutOff = 2
    #magCellSize = np.array([1,1,1])
    
    #generate all the parameters I will need
    h_list, k_list, l_list = generate_hkl(hkl_interval, direction)
    #test
    #h_list = 2*np.pi*h_list
    #k_list = 2*np.pi*k_list
    #l_list = 2*np.pi*l_list
    
    e_vals = generate_wt(omega_interval)
    
    mesh = np.zeros((len(h_list),len(e_vals)), dtype = complex)
    
    #mu_hat is a unit vector pointing in the direction of positive spin, which
    #I defined to be the direction of the spin of the atom whcih comes first in
    #the list.  I haven't though about whether the sign will make a difference
    #in the final result, but it does not really matter as the choice of
    #positive direction will always be arbitrary.
    spin = np.array(atom_list[0].spin)
    mu_hat = spin/np.sqrt(np.dot(spin,spin))
    
    ff_list = generate_form_factors(N_atoms_uc, atom_list, hkl_interval)
    s = []
    sigmas = []
    basis = []
    for a in atom_list:
        s.append(np.abs(np.dot(a.spin, mu_hat)))
        sigmas.append(np.dot(a.spin,mu_hat)/s[-1])
        #translate back to 0,0,0, but leave in crystallographic cell dimensions
        basis.append(a.pos - magCellSize)
        #I'm leaving it in crystallographic cell dimensions rather than magnetic
        #cell dimensions so that the periodicity of the dispersion relation
        #matches what we already calculate with the William's code, rather than
        #McQueeny's
        
    w_list = []
    Q_points = (np.array([h_list,k_list,l_list]).transpose())
    
    S_q = [] #For testing - S at the w points where S is maximized

    #calculate the cross section for each Q point (h,k,l)
    for Qindex in range(len(h_list)):
        Q = np.array(Q_points[Qindex])

        #create the matrix and find the eigenvalues and eigenvectors
        M = create_matrix(atom_list, jnums, jmats, numCutOff, magCellSize, Q-np.array(tau_list), mu_hat)
        w, evec = np.linalg.eig(M)
        
        w_list.append(w)#testing

        ff = ff_list[Qindex]

        #Get number of atoms
        numAtoms = numCutOff

        T = np.empty((numAtoms, numAtoms), dtype = complex)#normalized eigenvectors

        #Normalize evec to get T_ni
        for n in range(0, numAtoms):
            sgn_ln = w[n]/abs(w[n])
            tmp = 0.0
            for i in range(0, numAtoms):
                tmp += sigmas[i] * abs(evec[i][n])**2
            C = sgn_ln/tmp
            for i in range(0, numAtoms):
                T[i][n] = evec[i][n]*cmath.sqrt(C)
        #The indices of T actually go T[i][n], where n corresponds to the
        #eigenvalue, w[n], and i corresponds to the ith atom on the cell
        
        #1/12/12 The Matrix and T are tested for the simple antiferro case and
        #found to be correct
        
        #Now calculate the S(Q,w) / structure factor
        
        #mu_hat is a unit vector pointing in the direction of the spins
        #All the pins are aligned up to their sign, which does not matter here.
        #This method will work whether the spins are aligned along the z axis
        #or not.
        pol_factor = 0.5*(1.0 + (np.dot(mu_hat,Q)**2)/np.dot(Q,Q))
        
        gamma = ff #I need to make sure this is the same as his formfactor
        gamma = 1.0#testing
        
        #strfac = 0.0
        #For this one Q point, we get a whole S vs. e array
        tmp_e_dimension_array = np.zeros((len(e_vals)), dtype = complex)
        
        S_q.append(0.0)
        
        for n in range(numAtoms):
            strfac = 0.0
            for i in range(numAtoms):
                d_i = np.array(basis[i], dtype = np.float64)/np.array(magCellSize, dtype = np.float64)
                strfac += gamma*sigmas[i]*np.sqrt(s[i])*T[i][n]*np.exp(-1j*np.dot(Q,d_i))
            strfac = np.abs(strfac)**2
            #test:
            S_q[-1] += strfac
            #The 'Bose Occupation Factor,' as McQueeny calls it, which, I
            #believe is the Bose-Einstein distribution for identical bosons
            n_qn = 1.0/(np.exp(w[n]/(temp*kb)) - 1.0)
            
            #replace the delta function with a lorentzian
            for e_index in range(len(e_vals)):
                #delta(w - w_n(q))
                lorentzian1 = (1.0/np.pi)*HWHM/((e_vals[e_index] - w[n])**2 + HWHM**2)
                tmp_e_dimension_array += strfac*(n_qn+1.0)*lorentzian1
                #McQueeny left the other delta function out of his code, so for
                #now, I did too.
                
        #A test shows that numpy actually copies the values (not references)        
        mesh[Qindex] = pol_factor * tmp_e_dimension_array
                 
                 
    #Right now, mesh is only S(Q,w), not the full cross section, which is:
    #r_0^2 * k'/k * S(Q,w)
    r_0 = 1.0
    kp = 1.0
    k = 1.0             

    #print "Run Time: ", time.clock()-initial_time, " seconds"
    
    #testing
    w_list = np.array(w_list)
    for i in range(len(w_list)):
        print Q_points[i], "\t", S_q[i]
    plt.plot(np.linspace(hkl_interval[0],hkl_interval[1], hkl_interval[2])[2:-2], S_q[2:-2])
    plt.show()

    return (r_0**2) * (kp/k) * mesh

def get_Jij(jnums, jmats, num):
    """Returns the 3x3 numpy array corresponding to the matrix number, num from an atoms list
    of interactions."""
    index = jnums.index(num)
    return np.float64(np.array(jmats[index])) #This would be a sympy Matrix if I didn't use the np wrapper


def create_matrix(atom_list, jnums, jmats, numAtoms, magCellSize, q, mu_hat):
    """Using the Saenz technique.  The matrix indices, i and j refer to atoms in the magnetic
    unit cell, so M will be an N*N matrix, where N is the number of atoms in the magnetic cell.
    However, the l_vector is still in cyrstallographic unit cell units.

    -numCutOff is the number of atoms in the cutoff/interaction (or really magnetic, but we
    are assuming htey are the same) cell.
    -magCellSize is a length 3 np.array of the size of the magnetic cell in terms of crystallographic
    unit cells. ie. (1,1,1) for simple ferromagnet or (2,1,1) for simple antiferromagnetic chain along
    x direction.
    -q should also be a 3 element numpy array
    -atom_list, jums, and jmats come from readfiles.  atom_list must contain all atoms in the
    interaction/magnetic cell as well as atoms that they bond to.  Atoms in the magnetic cell
    should appear first in the list.
    -numAtoms is the number of atoms in the magnetic unit cell - this will be
    the size of the matrix.
    -mu_hat is a unit vector in the direction of the spins, which allows us to
    use any spin direction rather than the z axis
    
    This method works for the simple ferro and anti-ferro cases 1/12/12
    
    The units used for the vector, l, in this method are in terms of the
    magnetic unit cell, not the crystallographic.  I can't think of any sensible
    way to convert this to the crystallographic cell - the algorithm makes no
    sense then - it's not just affecting periodicity.  I think I may be able to
    simply divide out the unit cell at the end of the cross section calculation
    to change the periodicity, but not before the end.  If I use units of the
    magnetic cell here, I need to use those units for the rest of the calcuation
    so that they match.
    
    """
    magCellSize = np.array(magCellSize)#in case it's just a python array now

    M = np.zeros((numAtoms,numAtoms), dtype = np.complex)
    
    
    #Helper functions for creating the matrix
    def positions_match(atom1, atom2, magCellSize):
        pos1 = np.array(atom1.pos)
        pos2 = np.array(atom2.pos)
        #shift both positions back to first magCell - assuming they're positive
        while pos1[0] > magCellSize[0]:
            pos1[0] -= magCellSize[0]
        while pos1[1] > magCellSize[1]:
            pos1[1] -= magCellSize[1]
        while pos1[2] > magCellSize[2]:
            pos1[2] -= magCellSize[2]
        while pos2[0] > magCellSize[0]:
            pos2[0] -= magCellSize[0]
        while pos2[1] > magCellSize[1]:
            pos2[1] -= magCellSize[1]
        while pos2[2] > magCellSize[2]:
            pos2[2] -= magCellSize[2]
        pos_diff = np.array(pos2) - np.array(pos1)
        modulus = np.sqrt(np.dot(pos_diff, pos_diff))
        #This is a little sloppy, but I know I'm dealing with numbers on the
        #order of one, so I can use an absolute difference for comparison
        if modulus < 1e-8:
            return True
        return False
    
    def get_l(atom):
        #This is the 3-vector of integers indicating which magnetic cell the
        #atom is in.
        cell_size = np.array(magCellSize, dtype = np.int32)
        cell_0 = np.array(atom_list[0].pos, dtype = np.int32)/cell_size
        cell_l = np.array(atom.pos,dtype = np.int32)/cell_size
        return cell_l - cell_0
    
    


    #I'll do the first summation first - diagonal terms
    for i in range(numAtoms):
        for lk in range(len(atom_list[i].neighbors)):
            #assuming Jij matrix is diagonal, take one of the diagonal terms:
            J = (get_Jij(jnums, jmats,atom_list[i].interactions[lk]))[1][1]
            spin_k = np.dot(atom_list[atom_list[i].neighbors[lk]].spin,mu_hat)
            M[i][i] += J*spin_k
            
    #Now the second summation - which includes off-diagonal terms        
    for i in range(numAtoms):
        for j in range(numAtoms):
            spin_i = np.dot(atom_list[i].spin,mu_hat)
            spin_j = np.dot(atom_list[j].spin,mu_hat)
            S_i = np.abs(spin_i)
            S_j = np.abs(spin_j)
            sigma_j = spin_j/S_j
            #Now I want to find the jth atom in the lth cell
            #I'll lop through the interactions with i, and if any coincides
            #with the jth atom, I'll use that
            #A faster way would be to loop through the neighbors and find the
            #coinciding j for each.
            for nbr_index in range(len(atom_list[i].neighbors)):
                neighbor = atom_list[atom_list[i].neighbors[nbr_index]]
                if positions_match(atom_list[j], neighbor, magCellSize):
                    J = (get_Jij(jnums,jmats,atom_list[i].interactions[nbr_index]))[1][1]
                    M[i][j] -= sigma_j*np.sqrt(S_i*S_j)*J*np.exp(1j*np.dot(q,get_l(neighbor)))
    
    return M
        



#The following functions were copied from csection_calc, but also left in that
#file becuase they are used by both the MyQueeny algorithm and ours.  I
#duplicated the code becuase the code in this file is intended to be converted
#to Cython and I assume that talking between cython and python is slow.


#A copy of this was also left in cesction-calc because I assume talking between python and cyton is slow
def generate_hkl(hkl_interval, direction):

    #TODO - THIS DOES NOT WORK FOR ALL CASES! only when directional components are the same or 0
    #Fixed - testing

    #What if other points are 0 (not the first or last points)?
    #For that to happen they would need to generate an odd number of points across an interval -x to x
    if hkl_interval[0] == 0:
        hkl_interval[0] = 1e-6
    elif hkl_interval[1] == 0:
        hkl_interval[1] = -1e-6

    #print "hkl_interval: ", hkl_interval
    #print "direction: ", direction

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

#A copy of this was also left in cesction-calc because I assume talking between
#python and cython is slow
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
            #print el
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

#from csection_calc import LIFETIME_VALUE

if __name__ == "__main__":
    print "Running Test S(q)"
    #dir = '/home/tom/Desktop/spinwaves/'
    dir = '/media/0142-4351/Non-Ubuntu stuff/NIST/'
    csection(dir+'sg229-bcc-export.txt',dir+'sg229-bcc-spins.txt', [0,0,0], temp = 0.0001, HWHM = 1/np.pi,
            direction={'kx':1,'ky':1,'kz':1+1.0/8.0}, hkl_interval=[1e-12,12,1201],
            omega_interval=[0,5,100], kb = 0.086173423)
