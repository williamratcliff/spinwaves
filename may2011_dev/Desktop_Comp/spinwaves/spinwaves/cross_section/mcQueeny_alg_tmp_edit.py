import numpy as np
import time
from spinwaves.spinwavecalc.readfiles import atom, readFiles
from periodictable import elements


def mcQueeny_case(interactionfile, spinfile, tau_list, temperature,
                      direction={'kx':1,'ky':0,'kz':0}, hkl_interval=[1e-3,2*np.pi,100],
                      omega_interval=[0,5,100], eig_eps = 0.001):

    initial_time = time.clock()


    atom_list, jnums, jmats, N_atoms_uc, numCutOff, magCellSize = readFiles(interactionfile,spinfile, rtn_cutOffInfo = True)


    h_list, k_list, l_list = generate_hkl(hkl_interval, direction)
    e = generate_wt(omega_interval)
    ff_list = generate_form_factors(N_atoms_uc, atom_list, hkl_interval)
    s = []
    basis = []
    for a in atom_list:
        s.append(a.spin[2])#z component of spin (assume aligned with z)
        basis.append(a.pos - magCellSize) #translate  back to 0,0,0


    mesh = np.zeros((len(h_list), len(e)))

    for Qindex in range(len(h_list)):
        Q = (np.array([h_list,k_list,l_list]).transpose())[Qindex]

        M = create_matrix(atom_list, jnums, jmats, numCutOff, magCellSize, Q-np.array(tau_list))
        w, evec = np.linalg.eig(M)

        ff = ff_list[Qindex]


        numAtoms = numCutOff

        evec2 = evec.copy()

        T = evec
        sigma = []
        for i in range(numAtoms):
            sigma.append(s[i]/np.abs(s[i]))
        for n in range(numAtoms):
            norm = 0
            for i in range(numAtoms):
                norm += sigma[i]*(np.abs(evec2[n][i])**2)
            for i in range(numAtoms):
                T[n][i] = evec2[n][i]/np.sqrt(np.abs(norm))
        evec = T
        strfac = np.zeros(numAtoms)

        ff = 1

        for i in range(numAtoms):
            dum = 0
            for j in range(numAtoms):
                dum += (s[j]/np.sqrt(np.abs(s[j]))) * ff * evec[j][i] #* np.exp(1.0j*
                exp = np.exp(1.0j*(-Q[0]*basis[j][0] -Q[1]*basis[j][1] -Q[2]*basis[j][2]))
            strfac[i] = np.abs(dum)**2
        sqw = np.zeros(len(e))
        for j in range(numAtoms):
            for i in range(len(e)):
                y = LIFETIME_VALUE/2
                lorentzian = (1/np.pi)*y/((e[i]- w[j])**2+y**2)
                sqw[i] = sqw[i] + strfac[j] * lorentzian #he has some 'bose' factor in here
        mesh[Qindex] = sqw #* 0.5*(1 + Q[0]**2) #along x, Q_hat is just 1?

    print "Run Time: ", time.clock()-initial_time, " seconds"

def get_Jij(jnums, jmats, num):
    """Returns the 3x3 numpy array corresponding to the matrix number, num from an atoms list
    of interactions."""
    index = jnums.index(num)
    return np.float64(np.array(jmats[index])) #This would be a sympy Matrix if I didn't use the np wrapper

def generate_hkl(hkl_interval, direction):

    if hkl_interval[0] == 0:
        hkl_interval[0] = 1e-3
    elif hkl_interval[1] == 0:
        hkl_interval[1] = -1e-3

    print "hkl_interval: ", hkl_interval
    print "direction: ", direction
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


    return h_list, k_list, l_list

def generate_wt(omega_interval):

    return np.linspace(omega_interval[0],omega_interval[1],omega_interval[2])

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

        ff_list = Mq
    return np.array(ff_list)

from csection_calc import LIFETIME_VALUE