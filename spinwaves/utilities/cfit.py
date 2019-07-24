from __future__ import with_statement
import os
import sys
sys.path.append(os.path.abspath("../../../pak2010"))

import numpy as np
import sympy as sp

import park
from park.client import connect, JobDescription
from park.service.optimize import diffev, fitness

from cfit_utils import cfit_driver, numpy2sympy, read_files, read_tau
from spinwaves.cross_section.csection_calc import run_cross_section, run_eval_pointwise
from spinwaves.spinwavecalc.spinwave_calc_file import calculate_dispersion
from spinwaves.spinwavecalc.readfiles import atom, get_tokenized_line, findmat
from spinwaves.MonteCarlo.CSim import Sim_Aux, readFile, opt_aux


#    # Read interaction file
#    atoms, jMatrices = readFile(interaction_file)
#    jmats = matrices
#    k, tMax, tMin, tFactor = spins_data
#    # Calculate ground state spins from scratch
#    spins = Sim_Aux(k, tMax, tMin, tFactor, atoms, jmats)
#    jmats = numpy2sympy(jmats)
#    atomlist, N_atoms_uc = read_files(interaction_file)
#    tau_list = read_tau(tau_file)
#    # Pass in variable jMats and spins
#    # Get back 2D array of cross-section values
#    res = cfit_driver(jmats, atoms, spins, atomlist, N_atoms_uc, tau_list, 
#                      hkl_interval, w_interval, direction, tMin)
#    return res

@park.export
def fit_kernel(env, input):
    arr = np.load(input)
    items = (arr['atoms'], arr['spins'], arr['atomlist'], arr['N_atoms_uc'], 
             arr['tau_list'], arr['hkl_interval'], arr['w_interval'], 
             arr['direction'], arr['temperature'])
    (atoms, spins, atomlist, N_atoms_uc, tau_list, 
     hkl_interval, w_interval, direction, temperature) = items
#    jmats = lambda p: matrix_builder(p)
    return cfit_driver(lambda p: matrix_builder(p), atoms, spins, atomlist, N_atoms_uc, tau_list, 
         hkl_interval, w_interval, direction, temperature)

def matrix_builder(p):
    #get some list of parameters p = [J1, J2, J3, J4, J5, J6, J7, J8, J9]
    if len(p)%9:
        raise Exception("insufficient parameters to build matrices")
    matlist=[]
    i=0
    while i < len(p)/9:
        newmat = np.array([[p[i+0],p[i+1],p[i+2]],
                           [p[i+3],p[i+4],p[i+5]],
                           [p[i+6],p[i+7],p[i+8]]])
        matlist.append(newmat)
        i+=1
    return matlist

def fitting(matrices, interaction_file, tau_file, hkl_interval, w_interval, 
         direction, spins_data, save_file):
    """ 
    matrices must be list of sympy.matrices.Matrix 
    """
    atoms, jMatrices = readFile(interaction_file)
    jmats = matrices
    k, tMax, tMin, tFactor = spins_data
    spins = Sim_Aux(k, tMax, tMin, tFactor, atoms, jmats)
    jmats = numpy2sympy(jmats)
    atomlist, N_atoms_uc = read_files(interaction_file)
    tau_list = read_tau(tau_file)

    np.savez(save_file, atoms=atoms, spins=spins, atomlist=atomlist, 
             N_atoms_uc=N_atoms_uc, tau_list=tau_list, 
             hkl_interval=hkl_interval, w_interval=w_interval, 
             direction=direction, temperature=tMin)
    
    fit = dict(name="spinwaves.cross_section.fit_kernel", input=save_file, 
               files=dict(interactionfile="C:\\1_m.txt", 
                          spinfile="C:\\1_s.txt", 
                          taufile="C:\\1_t.txt", 
                          outfile="C:\\1_o.npz"))
    
    parameters = ('p8',0,-1,1), ('p7',0,-1,1), ('p6',0,-1,1), \
                 ('p5',0,-1,1), ('p4',0,-1,1), ('p3',0,-1,1), \
                 ('p2',0,-1,1), ('p1',0,-1,1), ('p0',0,-1,1)
    with connect("http://sparkle.ncnr.nist.gov:8000"):
        job = diffev(fit,parameters,ftol=1e-6,maxiter=100,npop=10)
    print job.wait()    

def main():
    save_file = "savefile.npz"
    interaction_file = "C:\\1_m.txt"
    tau_file = "C:\\1_t.txt"
    hkl_interval = [1e-3,2*np.pi,1000]
    w_interval = [0,5,1000]
    direction = {'kx':1,'ky':0,'kz':0}
    matrices = [np.array([[1,0,0],[0,0,0],[0,0,0]])]
    spins_data = [100, 10, 0.001, 0.9]#[k, tMax, tMin, tFactor]

    print matrices

    res = fitting(matrices, interaction_file, tau_file, hkl_interval, w_interval, 
               direction, spins_data, save_file)
    print res.shape

if __name__ == "__main__":
    main()
#    print matrix_builder([1,0,0,0,1,0,0,0,1,9])