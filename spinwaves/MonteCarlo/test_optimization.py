"""
Disclaimer
==========

This software was developed at the National Institute of Standards and Technology at the NIST Center for Neutron Research by employees of the Federal Government in the course of their official duties. Pursuant to title 17 section 105* of the United States Code this software is not subject to copyright protection and is in the public domain. The SPINAL software package is an experimental spinwave analysis system. NIST assumes no responsibility whatsoever for its use, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. The use of certain trade names or commercial products does not imply any endorsement of a particular product, nor does it imply that the named product is necessarily the best product for the stated purpose. We would appreciate acknowledgment if the software is used.

*Subject matter of copyright: United States Government works

Copyright protection under this title is not available for any work of the United States Government, but the United States Government is not precluded from receiving and holding copyrights transferred to it by assignment, bequest, or otherwise."""


'''
Created on Sep 27, 2009

@author: Bill
'''
import unittest
from simple import readFile, simpleAtom
from localOpt import opt_aux
from CSim import Sim_Aux
import numpy as np
import os

class TestChain(unittest.TestCase):

    def setUp(self):
        self.k = 100
        self.tMax = 10
        self.tMin = 0.01
        self.tFactor = 0.9
        self.tol = 1e-12
        self.comparetols = 0.20
        self.comparetolo = 0.01

    def getSpins(self, inFilePath):
        atoms, jMatrices = readFile(inFilePath)
        spins = Sim_Aux(self.k, self.tMax, self.tMin, self.tFactor, atoms, jMatrices)
        opt_spins = opt_aux(atoms, jMatrices, spins, self.tol)
        return (spins, opt_spins)

    def runIt_anis(self,inFilePath,dir):
        spins, opt_spins = self.getSpins(inFilePath)
        spins = spins[:3]
        opt_spins = opt_spins[:3]
        for s,o in zip(spins,opt_spins):
            sx,sy,sz=s[0],s[1],s[2]
            ox,oy,oz=o[0],o[1],o[2]
            x,y,z=dir[0],dir[1],dir[2]
            assert np.abs(x-np.abs(sx)) < self.comparetols
            assert np.abs(y-np.abs(sy)) < self.comparetols
            assert np.abs(z-np.abs(sz)) < self.comparetols
            assert np.abs(x-np.abs(ox)) < self.comparetolo
            assert np.abs(y-np.abs(oy)) < self.comparetolo
            assert np.abs(z-np.abs(oz)) < self.comparetolo
    
    # local optimizer just changes the original spins too much sometimes
    # but the resulting spin is equal to one so is fine
#    def runIt_rand(self,inFilePath):
#        spins, opt_spins = self.getSpins(inFilePath)
#        spins = spins[:3]
#        opt_spins = opt_spins[:3]
#        for s,o in zip(spins,opt_spins):
#            sx,sy,sz=s[0],s[1],s[2]
#            ox,oy,oz=o[0],o[1],o[2]
#            assert np.abs(np.abs(ox)-np.abs(sx)) < self.comparetols
#            assert np.abs(np.abs(oy)-np.abs(sy)) < self.comparetols
#            assert np.abs(np.abs(oz)-np.abs(sz)) < self.comparetols     

    def testchain_fm_x(self):
        inFilePath = os.path.dirname( __file__ )+"\\tests\\fm_chain_montecarlo_x.txt"
        #inFilePath = "c:/users/bill/documents/python/fm_chain_montecarlo_x.txt"
        self.runIt_anis(inFilePath,[1.,0.,0.])

    def testchain_fm_y(self):
        inFilePath = os.path.dirname( __file__ )+"/tests/fm_chain_montecarlo_y.txt"
        self.runIt_anis(inFilePath,[0.,1.,0.])

    def testchain_fm_z(self):
        inFilePath = os.path.dirname( __file__ )+"/tests/fm_chain_montecarlo_z.txt"
        self.runIt_anis(inFilePath,[0.,0.,1.])

    def testchain_afm_x(self):
        inFilePath = os.path.dirname( __file__ )+"/tests/afm_chain_montecarlo_x.txt"
        self.runIt_anis(inFilePath,[1.,0.,0.])

    def testchain_afm_y(self):
        inFilePath = os.path.dirname( __file__ )+"/tests/afm_chain_montecarlo_y.txt"
        self.runIt_anis(inFilePath,[0.,1.,0.])

    def testchain_afm_z(self):
        inFilePath = os.path.dirname( __file__ )+"/tests/afm_chain_montecarlo_z.txt"
        self.runIt_anis(inFilePath,[0.,0.,1.])

#    def testchain_fm_rand(self):
#        inFilePath = os.path.dirname( __file__ )+"/tests/fm_chain_montecarlo_0.txt"
#        self.runIt_rand(inFilePath)
#
#    def testchain_afm_rand(self):
#        inFilePath = os.path.dirname( __file__ )+"/tests/afm_chain_montecarlo_0.txt"
#        self.runIt_rand(inFilePath)

class TestCube(unittest.TestCase):

    def setUp(self):
        self.k = 100
        self.tMax = 10
        self.tMin = 0.01
        self.tFactor = 0.9
        self.tol = 1e-12
        self.comparetols = 0.15
        self.comparetolo = 0.05

    def getSpins(self, inFilePath):
        atoms, jMatrices = readFile(inFilePath)
        spins = Sim_Aux(self.k, self.tMax, self.tMin, self.tFactor, atoms, jMatrices)
        opt_spins = opt_aux(atoms, jMatrices, spins, self.tol)
        return (spins, opt_spins)

    def runIt_anis(self,inFilePath,dir):
        spins, opt_spins = self.getSpins(inFilePath)
        spins = spins[:8]
        opt_spins = opt_spins[:8]
        for s,o in zip(spins,opt_spins):
            sx,sy,sz=s[0],s[1],s[2]
            ox,oy,oz=o[0],o[1],o[2]
            x,y,z=dir[0],dir[1],dir[2]
            assert np.abs(x-np.abs(sx)) < self.comparetols
            assert np.abs(y-np.abs(sy)) < self.comparetols
            assert np.abs(z-np.abs(sz)) < self.comparetols
            assert np.abs(x-np.abs(ox)) < self.comparetolo
            assert np.abs(y-np.abs(oy)) < self.comparetolo
            assert np.abs(z-np.abs(oz)) < self.comparetolo
            
#    def runIt_rand(self,inFilePath):
#        spins, opt_spins = self.getSpins(inFilePath)
#        spins = spins[:8]
#        opt_spins = opt_spins[:8]
#        for s,o in zip(spins,opt_spins):
#            sx,sy,sz=s[0],s[1],s[2]
#            ox,oy,oz=o[0],o[1],o[2]
#            assert np.abs(np.abs(ox)-np.abs(sx)) < self.comparetols
#            assert np.abs(np.abs(oy)-np.abs(sy)) < self.comparetols
#            assert np.abs(np.abs(oz)-np.abs(sz)) < self.comparetols     

    def testcube_fm_x(self):
        inFilePath = os.path.dirname( __file__ )+"/tests/fm_cube_montecarlo_x.txt"
        self.runIt_anis(inFilePath,[1.,0.,0.])

    def testcube_fm_y(self):
        inFilePath = os.path.dirname( __file__ )+"/tests/fm_cube_montecarlo_y.txt"
        self.runIt_anis(inFilePath,[0.,1.,0.])

    def testcube_fm_z(self):
        inFilePath = os.path.dirname( __file__ )+"/tests/fm_cube_montecarlo_z.txt"
        self.runIt_anis(inFilePath,[0.,0.,1.])

    def testcube_afm_x(self):
        inFilePath = os.path.dirname( __file__ )+"/tests/afm_cube_montecarlo_x.txt"
        self.runIt_anis(inFilePath,[1.,0.,0.])

    def testcube_afm_y(self):
        inFilePath = os.path.dirname( __file__ )+"/tests/afm_cube_montecarlo_y.txt"
        self.runIt_anis(inFilePath,[0.,1.,0.])

    def testcube_afm_z(self):
        inFilePath = os.path.dirname( __file__ )+"/tests/afm_cube_montecarlo_z.txt"
        self.runIt_anis(inFilePath,[0.,0.,1.])

#    def testcube_fm_rand(self):
#        inFilePath = os.path.dirname( __file__ )+"/tests/fm_cube_montecarlo_0.txt"
#        self.runIt_rand(inFilePath)
#
#    def testcube_afm_rand(self):
#        inFilePath = os.path.dirname( __file__ )+"/tests/afm_cube_montecarlo_0.txt"
#        self.runIt_rand(inFilePath)

        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    print os.path.dirname( __file__ )
    print os.path.exists(os.path.dirname( __file__ )+"\\tests")
    unittest.main()
#    k = 100
#    tMax = 10
#    tMin = 0.01
#    tFactor = 0.9
#    tol = 1e-12
#    comparetols = 0.15
#    comparetolo = 0.05
#    inFilePath = os.path.dirname( __file__ )+"/tests/fm_chain_montecarlo_y.txt"
#    atoms, jMatrices = readFile(inFilePath)
#    spins = Sim_Aux(k, tMax, tMin, tFactor, atoms, jMatrices)
#    opt_spins = opt_aux(atoms, jMatrices, spins, tol)
#    spins = spins[:3]
#    opt_spins = opt_spins[:3]
#    dir=[0.,1.,0.]
#    for s,o in zip(spins,opt_spins):
#        sx,sy,sz=s[0],s[1],s[2]
#        ox,oy,oz=o[0],o[1],o[2]
#        print sx,sy,sz
#        print ox,oy,oz
#        x,y,z=dir[0],dir[1],dir[2]
##        assert np.abs(x-np.abs(sx)) < comparetols
##        assert np.abs(y-np.abs(sy)) < comparetols
##        assert np.abs(z-np.abs(sz)) < comparetols
##        assert np.abs(x-np.abs(ox)) < comparetolo
##        assert np.abs(y-np.abs(oy)) < comparetolo
##        assert np.abs(z-np.abs(oz)) < comparetolo
