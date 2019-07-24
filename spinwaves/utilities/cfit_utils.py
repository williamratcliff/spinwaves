"""
The purpose of this module is to create a single function which will calculate
the cross-section starting from the simulated annealing process. The problem 
with the current routine is that many files are saved, then read in that process
and it is difficult to write a fitter around the routine.
"""

"""
input needs to contain:
-- k, tMax, tMin, tFactor, atoms, jMatrice for CSIM
-- interaction, spin files for SWC
-- interaction, spin, tau, out files for CS
-- hkl range/steps, omega range/steps, temp for CSIM

interaction, spin files -> atoms
k - csim steps
tMax, tMin - csim temp range
tFactor - csim temp factor
"""
# import
import os
import sys
import time

import sympy as sp
import numpy as np

from spinwaves.cross_section.csection_calc import run_cross_section, run_eval_pointwise
from spinwaves.spinwavecalc.spinwave_calc_file import calculate_dispersion
from spinwaves.spinwavecalc.readfiles import atom, get_tokenized_line, findmat
from spinwaves.MonteCarlo.CSim import Sim_Aux, readFile, opt_aux

def read_tau(tau_file):

    failed = False
    tau_text = ''
    tau_list = []

    try:
        tau_file = open(tau_str,'r')
        tau_text = tau_file.read()
    except:
        failed = True

    items = tau_text.split()
    if len(items)%3 and not len(items):
        failed = True

    i = 0
    while not failed and i <= len(items)-3:
        tau1, tau2, tau3 = None, None, None
        try:
            tau1 = float(items[i])
            tau2 = float(items[i+1])
            tau3 = float(items[i+2])
        except:
            failed = True
        tau_list.append([tau1,tau2,tau3])
        i+=3
    return tau_list

def add_spins(atoms, spins, atomlist):
    for atom1 in atomlist:
        for i in range(len(atoms)):
            xPos = atoms[i].pos[0]
            yPos = atoms[i].pos[1]
            zPos = atoms[i].pos[2]
            if (atom1.pos[0] == xPos and atom1.pos[1] == yPos and atom1.pos[2] == zPos):
                spin = spins[i]
                atom1.spin = spin
                rmat = findmat(spin)
                atom1.spinRmatrix = rmat
                break
    return atomlist

def read_files(interactionFileStr, allAtoms=False):
    """ MODIFIED FROM SPINWAVE CALC FILE TO ONLY READ INTERACTION FILE"""
    """modified from read_interactions.  Originally this(read_interactions) read
    in the atoms from the interaction file and matched the spin rotation
    matrices with the appropriate atoms based on indices.  Now it takes the
    interaction and spin file strings, reads in the atom information and matches
    it with spin rotation matrices based on coordinates"""
    #print interactionFileStr
    interactionFile = open(interactionFileStr, 'r')
    myFlag=True
    returnline=['']
    jmats=[]
    jnums=[]
    atomlist=[]
    numcell=0
    Na, Nb, Nc = 0,0,0
    while myFlag:
        tokenized=get_tokenized_line(interactionFile,returnline=returnline)
        if not(tokenized):
            break
        if tokenized==[]:
            break
        
        if tokenized[0] == "#na":
            tokenized=get_tokenized_line(interactionFile,returnline=returnline)
            Na = int(tokenized[0])
            Nb = int(tokenized[1])
            Nc = int(tokenized[2])
            
        if tokenized[0]=='#number':
            while 1:
                tokenized=get_tokenized_line(interactionFile,returnline=returnline)
                #print 'intoken ',tokenized
                if tokenized==[]:
                    break
                if tokenized[0]!='#atomnumber':
                    #print tokenized[0]
                    jnum=float(tokenized[0])
                    j11=float(tokenized[1])
                    j12=float(tokenized[2])
                    j13=float(tokenized[3])
                    j21=float(tokenized[4])
                    j22=float(tokenized[5])
                    j23=float(tokenized[6])
                    j31=float(tokenized[7])
                    j32=float(tokenized[8])
                    j33=float(tokenized[9]) 
                    #jij=N.matrix([[j11,j12,j13],[j21,j22,j23],[j31,j32,j33]],'Float64')
                    jij=sp.matrices.Matrix([[j11,j12,j13],[j21,j22,j23],[j31,j32,j33]])
                    jnums.append(jnum)
                    jmats.append(jij)
                else:
                    currnum=0
                    while 1:
                        tokenized=get_tokenized_line(interactionFile,returnline=returnline)
                        if not(tokenized):
                            break
                        #print tokenized
                        atom_num=int(tokenized[0])
                        if tokenized[5] == 'x' or allAtoms:  #If it is in the first interaction cell
                            print "atom in first interaction cell"
                            x,y,z=float(tokenized[6]),float(tokenized[7]),float(tokenized[8])
                            print x,y,z

                            Dx,Dy,Dz=float(tokenized[9]),float(tokenized[10]),float(tokenized[11])
                            #spin0=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'float64')
                            spinMagnitude = float(tokenized[12])
                            valence = int(tokenized[4])
                            massNum = int(tokenized[3])
                            label = tokenized[1].capitalize()
                            pos0=[x,y,z]
                            atom0=atom(spinMag=spinMagnitude, pos=pos0,label = label, Dx=Dx,Dy=Dy,Dz=Dz, orig_Index = atom_num, valence = valence, massNum = massNum)
                            neighbors=[]
                            interactions=[]
                            for i in range(13,len(tokenized),2):
                                interacting_spin=int(tokenized[i])
                                #index number in export list not necessarily the same as index
                                #in list of atoms in first interacting 'cell'
                                interaction_matrix=int(tokenized[i+1])
                                neighbors.append(interacting_spin)
                                interactions.append(interaction_matrix)

                            atom0.neighbors=neighbors
                            atom0.interactions=interactions
                            currnum=currnum+1
                            atomlist.append(atom0)
    interactionFile.close()
    
    def inDesiredCell(atom):
        if atom.pos[0] >= Na and atom.pos[0] < (Na + 1):
            if atom.pos[1] >= Nb and atom.pos[1] < (Nb + 1):
                if atom.pos[2] >= Nc and atom.pos[2] < (Nc + 1):
                    return True
        return False
        
    #Add atoms in desired cell to the beginning of the new list
    newAtomList = []
    Flag=True
    i=0
    while Flag:
    #for i in range(len(atomlist)):
        if inDesiredCell(atomlist[i]):
            numcell += 1
            newAtomList.append(atomlist.pop(i))
        else:
            i=i+1
        if i==len(atomlist):
            Flag=False
    
    #Add remaining atoms to the new list
    for i in range(len(atomlist)):
        newAtomList.append(atomlist[i])
    atomlist = newAtomList

    for a in atomlist:
        neighborList = a.neighbors
        i = 0
        while i < len(neighborList):
            neighbor_index = neighborList[i]
            for atom_index in range(len(atomlist)):
                atom1 = atomlist[atom_index]
                if atom1.origIndex == neighbor_index:
                    a.neighbors[i] = atom_index
                    i +=1
                    break
            else:
                neighborList.pop(i)
                a.interactions.pop(i)
            
    for atom1 in atomlist:
        print atom1.pos, atom1.neighbors
    
    return atomlist, numcell

def cfit_driver(jmats, atoms, spins, atomlist, N_atoms_uc, tau_list, 
         hkl_interval, w_interval, direction, temperature):
    
    # Ground State
    opt_spins = opt_aux(atoms, jmats, spins)
    atom_list = add_spins(atoms, spins, atomlist)
    N_atoms = len(atom_list)
    
    # Dispersion
    Hsave = calculate_dispersion(atom_list, N_atoms_uc, N_atoms, jmats, 
                                 showEigs = False)
    
    # Cross-section
    (atom_list, N_atoms_uc, csection, 
     kaprange, tau_list, eig_list, 
     kapvect, wt_list, fflist) = run_cross_section(None, None, tau_list, temperature,
                                                   direction, hkl_interval, w_interval, 
                                                   atom_list=atom_list, N_atoms_uc=N_atoms_uc,
                                                   Hsave=Hsave, nosave=True)
    x,y,z = run_eval_pointwise(N_atoms_uc, atom_list, csection, kaprange, 
                               tau_list, eig_list, kapvect, wt_list, fflist, 
                               temperature, direction)
    
    # Return 2d array
    return z

def numpy2sympy(matrices):
    mats = []
    for mat in matrices:
        newmat = sp.matrices.Matrix([[mat[0][0],mat[0][1],mat[0][2]],
                                     [mat[1][0],mat[1][1],mat[1][2]],
                                     [mat[2][0],mat[2][1],mat[2][2]]])
        mats.append(newmat)
    return mats

def main(matrices, interaction_file, tau_file, hkl_interval, w_interval, 
         direction, spins_data):
    """ 
    matrices must be list of sympy.matrices.Matrix 
    """
    # Read interaction file
    atoms, jMatrices = readFile(interaction_file)
    jmats = matrices
    k, tMax, tMin, tFactor = spins_data
    # Calculate ground state spins from scratch
    spins = Sim_Aux(k, tMax, tMin, tFactor, atoms, jmats)
    jmats = numpy2sympy(jmats)
    atomlist, N_atoms_uc = read_files(interaction_file)
    tau_list = read_tau(tau_file)
    # Pass in variable jMats and spins
    # Get back 2D array of cross-section values
    res = cfit_driver(jmats, atoms, spins, atomlist, N_atoms_uc, tau_list, 
                      hkl_interval, w_interval, direction, tMin)
    return res

if __name__ == "__main__":
    # Fake data
    interaction_file = "C:\\1_m.txt"
    tau_file = "C:\\1_t.txt"
    hkl_interval = [1e-3,2*np.pi,1000]
    w_interval = [0,5,1000]
    direction = {'kx':1,'ky':0,'kz':0}
    matrices = [np.array([[1,0,0],[0,0,0],[0,0,0]])]
    spins_data = [100, 10, 0.001, 0.9]#[k, tMax, tMin, tFactor]

    print matrices

    res = main(matrices, interaction_file, tau_file, hkl_interval, w_interval, 
               direction, spins_data)
    print res.shape
    
