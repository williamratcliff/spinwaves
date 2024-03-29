"""
Disclaimer
==========

This software was developed at the National Institute of Standards and Technology at the NIST Center for Neutron Research by employees of the Federal Government in the course of their official duties. Pursuant to title 17 section 105* of the United States Code this software is not subject to copyright protection and is in the public domain. The SPINAL software package is an experimental spinwave analysis system. NIST assumes no responsibility whatsoever for its use, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. The use of certain trade names or commercial products does not imply any endorsement of a particular product, nor does it imply that the named product is necessarily the best product for the stated purpose. We would appreciate acknowledgment if the software is used.

*Subject matter of copyright: United States Government works

Copyright protection under this title is not available for any work of the United States Government, but the United States Government is not precluded from receiving and holding copyrights transferred to it by assignment, bequest, or otherwise."""


import string
import copy
import xml.dom.minidom
#import xml.dom.ext
import wx.grid
from wx.py.dispatcher import send, connect
import spinwaves.vtkModel.CifFile
#import vtkModel.CifFile as CifFile
from spinwaves.vtkModel import SpaceGroups
from spinwaves.vtkModel.CellClass import Cell
from spinwaves.vtkModel.MagneticCellClass import MagneticCell
import numpy
import time
import datetime
from spinwaves.vtkModel.BondClass import JParam
from spinwaves.vtkModel.Parameter_Manager import ParamManager
import spinwaves.vtkModel.CifFile as CifFile
from spinwaves.vtkModel.SpaceGroups import GetSpaceGroup
from spinwaves.vtkModel.AtomClass import getElObj
from spinwaves.vtkModel.SymmetryUtilities import equalPositions
from spinwaves.vtkModel.AtomClass import getElObj

class Session():
    """Stores information about a user session

    This class stores and manipulates the model (magnetic or Cuttof Cell).The
    GUI stores the instance of Session that is in use, but to communicate with
    the GUI, Session sends messages.(The GUI is aware of Session, but Session
    is not aware of the GUI.)"""
    def __init__(self):
        self.parameter_manager = ParamManager()
        #Stores bond information
        self.bondTable = bondTable(self.parameter_manager)
        #Stores atom information
        self.atomTable = atomTable()
        self.MagCell = None
        #For now the magnetic cell class will be used for the cutoff cell
        #receive message from the fitresultwindow to copy the fit parameters in the bondTable
        connect(self.OnUseFitData, signal = "Use Fit Data")
        
    def OnUseFitData(self, bondTable):
        self.bondTable.data = bondTable.data
        print self.bondTable.data[0][8][0][0].value
        send(signal = "Model Change", sender = "Session")
        
    def getAtomTable(self):
        return self.atomTable
    
    def getBondTable(self):
        return self.bondTable
    
    def getMagneticCell(self):
        return self.MagCell
    
    #For now magnetic cell and cuttoff cell are exactly the same so magnetic cell will
    #be treated as cutoff cell
    def getCutoffCell(self):
        return self.MagCell
    
    def loadXMLStr(self, str, notifyGUI = False):
        """this loads all the information in an xml string to the model."""
        doc = xml.dom.minidom.parseString(str)
    
        #Get Main node
        nodes = doc.getElementsByTagName('Spinwaves_Session')
        if len(nodes) == 1:
            mainNode = nodes[0]
        else:#There should only be one of these elements
            raise Exception('not valid XML Session file')
        
        #Get unit cell node
        nodes = mainNode.getElementsByTagName('unit_cell')
        if len(nodes) == 1:
            unitCellNode = nodes[0]
        else:#There should only be one of these elements
            raise Exception('not valid XML Session file')
        
        #Get magnetic cell node
        nodes = mainNode.getElementsByTagName('magnetic_cell')
        if len(nodes) == 1:
            magCellNode = nodes[0]
        else:#There should only be one of these elements
            raise Exception('not valid XML Session file')
        
        #Get atoms node
        nodes = mainNode.getElementsByTagName('atoms')
        if len(nodes) == 1:
            atomsNode = nodes[0]
        else:#There should only be one of these elements
            raise Exception('not valid XML Session file')
        
        #Get bonds node
        nodes = mainNode.getElementsByTagName('bonds')
        if len(nodes) == 1:
            bondsNode = nodes[0]
        else:#There should only be one of these elements
            raise Exception('not valid XML Session file')
        
        
        spaceGroupInt = int(unitCellNode.getAttribute('space_group'))
        a = float(unitCellNode.getAttribute('a'))
        b = float(unitCellNode.getAttribute('b'))
        c = float(unitCellNode.getAttribute('c'))
        alpha = float(unitCellNode.getAttribute('alpha'))
        beta = float(unitCellNode.getAttribute('beta'))
        gamma = float(unitCellNode.getAttribute('gamma'))
        Na = int(magCellNode.getAttribute('Na'))
        Nb = int(magCellNode.getAttribute('Nb'))
        Nc = int(magCellNode.getAttribute('Nc'))
        
        #Get atoms
        atomNodes = atomsNode.getElementsByTagName('atom')
        atomData = []
        for i in range(len(atomNodes)):
            name = atomNodes[i].getAttribute('Name')
            #I called this massNum, but it's really the atomic number now
            massNum = atomNodes[i].getAttribute('Atomic_Number')
            x = atomNodes[i].getAttribute('x')
            y = atomNodes[i].getAttribute('y')
            z = atomNodes[i].getAttribute('z')
            Dx = atomNodes[i].getAttribute('SIA_Dx')
            Dy = atomNodes[i].getAttribute('SIA_Dy')
            Dz = atomNodes[i].getAttribute('SIA_Dz')
            spinMag = atomNodes[i].getAttribute('spin_magnitude')
            valence = atomNodes[i].getAttribute('valence')
            
            atomData.append([name, int(massNum), float(x),float(y),float(z),
                             float(Dx), float(Dy), float(Dz), float(spinMag),
                             valence])
            
            self.atomTable.SetValue(i, 0, name)
            self.atomTable.SetValue(i, 1, massNum)
            self.atomTable.SetValue(i, 2, valence)
            self.atomTable.SetValue(i, 3, x)
            self.atomTable.SetValue(i, 4, y)
            self.atomTable.SetValue(i, 5, z)
            self.atomTable.SetValue(i, 6, Dx)
            self.atomTable.SetValue(i, 7, Dy)
            self.atomTable.SetValue(i, 8, Dz)
            self.atomTable.SetValue(i, 9, spinMag)

        
        bondNodes = bondsNode.getElementsByTagName('bond')
#        bondData = []
        
        #parameters are split up by the J matrix they are in, but need to be resolved
        #into a single list.  This method is relying on the parameter numbers starting
        #at 1 and not having any gaps, which should always be the case
        allParamNodes = []
        
        for i in range(len(bondNodes)):
            atom1Num = bondNodes[i].getAttribute('Atom1_Number')
            atom1Na = bondNodes[i].getAttribute('Atom1_Na')
            atom1Nb = bondNodes[i].getAttribute('Atom1_Nb')
            atom1Nc = bondNodes[i].getAttribute('Atom1_Nc')
            atom2Num = bondNodes[i].getAttribute('Atom2_Number')
            atom2Na = bondNodes[i].getAttribute('Atom2_Na')
            atom2Nb = bondNodes[i].getAttribute('Atom2_Nb')
            atom2Nc = bondNodes[i].getAttribute('Atom2_Nc')
            on = bondNodes[i].getAttribute('On')
            self.bondTable.SetValue(i, 0, atom1Num)
            self.bondTable.SetValue(i, 1, atom1Na)
            self.bondTable.SetValue(i, 2, atom1Nb)
            self.bondTable.SetValue(i, 3, atom1Nc)
            self.bondTable.SetValue(i, 4, atom2Num)
            self.bondTable.SetValue(i, 5, atom2Na)
            self.bondTable.SetValue(i, 6, atom2Nb)
            self.bondTable.SetValue(i, 7, atom2Nc)
            self.bondTable.SetValue(i, 9, on)
            
            #Jmatrix 
            nodes = bondNodes[i].getElementsByTagName('jMatrix')
            #if len(nodes) > 0:   
            if len(nodes) == 1:
                jNode = nodes[0]
            else:#There should only be one of these elements
                raise Exception('not valid XML Session file')
            
            parameterNodes = jNode.getElementsByTagName('Parameter')
            if parameterNodes.length!=9:
                raise Exception('not valid XML Session file: There should be 9 parameters per J-matrix')
            

            #Add the nodes to the full list of parameters
            for paramNode in parameterNodes:
                allParamNodes.append(paramNode)
            
        #insert the parameters into the manager in the correct order
        manager = ParamManager()
        tiedLists = []
        for i in range(len(allParamNodes)):
            #i is also the current index in the manager list of parameters
            for node in allParamNodes:
                if int(node.getAttribute('name')[1:]) == i:
                    fitStr = node.getAttribute('fit')
                    if fitStr =='True':
                        fit = True
                    else:
                        fit = False
                    value = float(node.getAttribute('value'))
                    min = node.getAttribute('min')
                    max = node.getAttribute('max')
                    default = float(node.getAttribute('default'))
                    
                    param = JParam(manager, fit, value, min , max, default)
                    print len(manager.parameters)
                    tiedToStr = node.getAttribute('tiedTo')
                    strList = tiedToStr[1:len(tiedToStr)-1].split(',')
                    tied = []
                    for tie in strList:
                        if tie != '':
                            tied.append(int(tie))
                    tiedLists.append(tied)
                    break
            else:
                raise Exception("Incorrect parameter indices in XML file.")
            
        for i in range(len(manager.parameters)):
            manager.parameters[i].tieToMany(tiedLists[i])
        self.parameter_manager = manager
        self.bondTable.paramManager = manager
        
        #now construct the matrices with the appropriate JParam objects
        numMatrices = len(allParamNodes)/9
        print 'matrice: ', numMatrices
        def findParam(posName, startIndex, endIndex):
            """find the jParam object in the node list between indices startindex and endIndex.
            (endIndex not included)"""
            for i in range(startIndex, endIndex):
                if allParamNodes[i].getAttribute('position') == posName:
                    return manager.parameters[i]
            
            raise Exception('Problem resolving parameter to J Matrix from xml file')
        
        for i in range(numMatrices):
            j11 = findParam('j11', i*9, i*9+9)
            j12 = findParam('j12', i*9, i*9+9)
            j13 = findParam('j13', i*9, i*9+9)
            j21 = findParam('j21', i*9, i*9+9)
            j22 = findParam('j22', i*9, i*9+9)
            j23 = findParam('j23', i*9, i*9+9)
            j31 = findParam('j31', i*9, i*9+9)
            j32 = findParam('j32', i*9, i*9+9)
            j33 = findParam('j33', i*9, i*9+9)
            mat = numpy.array([[j11, j12, j13],
                               [j21, j22, j23],
                               [j31, j32, j33]])
            
            self.bondTable.SetValue(i, 8, mat)

        
        #self.cellChange(spaceGroupInt, a, b, c, alpha, beta, gamma, Na, Nb, Nc, Na, Nb, Nc, atomData, notifyGUI)
        self.updateCell(spaceGroupInt, a, b, c, alpha, beta, gamma, Na, Nb, Nc,
                        Na, Nb, Nc, atomData)
        if notifyGUI:
            self.refreshGUI()
        
        return spaceGroupInt, a, b, c, alpha, beta, gamma, Na, Nb, Nc, Na, Nb, Nc, atomData
        
    def openXMLSession(self, filename):
        """This loads an xml session file by calling loadXMLStr and then notifying the GUI."""
        #doc = xml.dom.minidom.parse(filename)
        handle = open(filename)
        xmlStr = handle.read()
        
        spaceGroupInt, a, b, c, alpha, beta, gamma, Na, Nb, Nc, Na, Nb, Nc, atomData = self.loadXMLStr(xmlStr, notifyGUI = True)
        
        #Send Message to GUI
        send(signal = "File Load", sender = "Session", spaceGroup = spaceGroupInt, a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma, magNa = Na, magNb = Nb, magNc = Nc, cutNa = Na, cutNb = Nb, cutNc = Nc)
        send(signal = "Generate Bonds", sender = "Session", event=None)
        
    def clearBonds(self):
        self.MagCell.clearAllBonds()
                
    #def changeBonds(self, bondData, notifyGUI = True):
        #"""This method is called when the bonds need to be changed in the model,
        #such as when info is entered in the GUI or when a file is loaded.
        #If bondData includes two bonds connecting the same two atoms, an
        #exception will be raised."""
        ##This function was implemented before addBond() or refreshGUI()
        #self.MagCell.clearAllBonds()
        #for i in range(len(bondData)):
            #if len(bondData[i]) < 10 or bondData[i][9] == 'X': #When this method is called by the load file method, it includes the rows that are off
                #cell1 = self.MagCell.cellAtPosition((int(bondData[i][1]), int(bondData[i][2]), int(bondData[i][3])))
                #cell2 = self.MagCell.cellAtPosition((int(bondData[i][5]), int(bondData[i][6]), int(bondData[i][7])))
                
                #atom1 = cell1.atomAtIndex(int(bondData[i][0]) - 1)
                #atom2 = cell2.atomAtIndex(int(bondData[i][4]) - 1)
                
    
                ##bondData[i][8] can be '' instead of None
                #if not isinstance(bondData[i][8], str):
                    #if not self.MagCell.hasBond(atom1, atom2, bondData[i][8]):
                        #self.MagCell.addBond(atom1, atom2, bondData[i][8])
                #else:
                    #if not self.MagCell.hasBond(atom1, atom2):
                        #self.MagCell.addBond(atom1, atom2)
                        
        #if notifyGUI:
            #send(signal = "Model Change", sender = "Session")
                        
    def addBond(self, **kwds):
        """This will add a bond connecting atom1 and atom2.
	
        Keywords:	
	(atom1Num, atom1CellPos, atom2Num, atom2CellPos, jMatrix, atom1, atom2)
        
        Either the atom or it's number and cell position must be supplied.
	
	atom1Num is the unique number identifying the given atom in the unit cell.
	atom1CellPos is a 3 integer tuple with the coordinates of the unit
	cell containing atom1.
        
        JijMatrix is the 3x3 matrix describing the spin-spin interaction between
        the atoms for hte heisenberg hamiltonian.
        
        If a bond with the same JMatrix that links the same two atoms already exists 
        nothing happens.  If a bond linking the two atoms exists, but has a different
        JMatrix, an Exception is raised.
	"""
        print kwds
        self.MagCell.addBond(**kwds)
        
        
    def newCell(self,spaceGroupInt,a,b,c,alpha, beta, gamma, magNa, magNb, magNc, cutNa, cutNb, cutNc):
        """This method is called when a change needs to be made to the
        crystallographic or cutoff(or magnetic), such as when the user enters
        info in the GUI or when a file is loaded.  To add atoms to the new
        unit cell, call addAtom().
        """
        spaceGroup = SpaceGroups.GetSpaceGroup(spaceGroupInt)
        unitcell = Cell(spaceGroup, 0,0,0, a, b, c, alpha, gamma, beta)
        self.MagCell = MagneticCell(unitcell, cutNa, cutNb, cutNc, spaceGroup)
        
        #Regenerate Bonds as well
        #self.changeBonds(self.bondTable.data, notifyGUI)
        
    def addAtom(self, symbol, position, massNum = None, **kwds):
        """Given the information of one atom and the space group associated with
        this Cell object, this method creates all the symmetry equivalent atoms
        and adds them to the model (cutoff cell).
        
        -symbol is the element symbol (H, He, Li...)
	-massNum is the atomic Mass Number for the element
	If the symbol and mass number do not correspond to a real element, an
	ElementNotFound Exception will be raised.
	
        -position is the fractional coordinates in the unit cell (a,b,c)
        
        Optional keyword arguments are:
        (spin, spinMagnitude, valence, anisotropy, rgb, radius, and description)
        
        -spin is a tuple(Sx, Sy, Sz). It is optional.
        -spinMagnitude is the total magnitude of the spin =sqrt(Sx^2+Sy^2+Sz^2)
        -anisotropy is the single ion anisotropy of the atom (Dx, Dy, Dz)
        -rgb is a 3 element tuple or list describing the color of the atom.  If
        this is notsupplied, the default color for that element is used.
        -radius is the radius of the atom in angstroms.  If it is None, it will
        be given the default value for the given element.
        -description is a string describing the atom such as a name/label and is
        optional.
        -Valence is a string describing the charge of the atom.
        """
        self.MagCell.addAtom(symbol, position, massNum, **kwds)
        
    def updateCell(self, spaceGroupInt, a,b,c, alpha, beta, gamma, magNa, magNb,
                   magNc, cutNa, cutNb, cutNc, atomData):
        """Updtates the cell with the supplied parameters.
        
        The cell is recreated and the atoms in atomData recreated.  If this is
        the first time any cutoff cell has been created, then a newone will be
        made.  Otherwise, bonds that can be preserved will be if the atoms that
        it connects still exist.  This will be determined based on position.
        If the two ID numbers associated with the atoms of the bond refer to 
        atoms whose positions have not changed, then the bond will also be
        recreated.  Any bonds whose ID numbers refer to nonexistant atoms
        or atoms whose positions have changed will be disabled.  The bond
        information in the bondtable will not be changed.
        
        atomData is a list where each entry corresponds to one atom and is a
        list of it's parameters in the following format:
        [Symbol,MassNumber,aPos,bPos,cPos,Da,Db,Dc,spinMagnitude,valence]
                                           ^anisotropy
        """
        oldBonds = None
        if self.MagCell != None:
            oldBonds = self.MagCell.getBonds()
            
        self.newCell(spaceGroupInt,a,b,c,alpha, beta, gamma, magNa, magNb,
                     magNc, cutNa, cutNb, cutNc)
        for data in atomData:
            self.addAtom(symbol = data[0], massNum = data[1],#massNum will most likely be None
                         position = (data[2], data[3], data[4]),
                         anisotropy = (data[5], data[6], data[7]),
                         spinMagnitude = data[8], valence = data[9])
    
        
        if oldBonds != None:  #There were bonds that need be sorted out
            for oldBond in oldBonds:
                atom1 = oldBond.getAtom1()
                atom2 = oldBond.getAtom2()
                pos1 = atom1.getPosition()
                pos2 = atom2.getPosition()
                cell1 = self.MagCell.cellAtPosition(pos1)
                cell2 = self.MagCell.cellAtPosition(pos2)
                if cell1 != None and cell2 != None:
                    newAtom1 = cell1.atomWithID(atom1.getIDNum())
                    newAtom2 = cell2.atomWithID(atom2.getIDNum())
                    if newAtom1 != None and newAtom2 != None:
                        if equalPositions(newAtom1.getPosition(), pos1) and \
                           equalPositions(newAtom2.getPosition(), pos2):
                            #Save this Bond
                            #Error condition should be impossible
                            #Becuase many of these bonds are likely symmetry
                            #equivalent, this method will be called for bonds
                            #that already exist, but they won't be created.
                            id1 = newAtom1.getIDNum()
                            id2 = newAtom2.getIDNum()
                            self.addBond(atom1Num = id1,
                                         atom1CellPos = cell1.getPosition(),
                                         atom2Num = id2,
                                         atom2CellPos = cell2.getPosition())
        
        #Nothing is done to the bondtable for now   
        #self.refreshGUI()
        
        
    def refreshGUI(self):
        send(signal = "Model Change", sender = "Session")
        
    #defcellChange(self,spaceGroupInt,a,b,c,alpha, beta, gamma, magNa, magNb, magNc, cutNa, cutNb, cutNc, atomData, notifyGUI = True):
        #"""This method is called when a change needs to be made to the
        #crystallographic or cutoff(or magnetic) cells or the atoms in them, such
        #as when the user enters info in the GUI or when a file is loaded.

        #atomData format (each line):
        #label massNum aPos bPos cPos anisotropy_a anisotropy_b anistropy_c spinMag valence spinTuple elementObj
                
        #elementObj is optional.  If each entry in the list has a length of 11(elementObj not included)
        #then label must be the symbol of the element(H, He, Li...).  If elementObj is not included
        #and label does not correspond to an element symbol, and Exception will be generated.
        
        #spinTuple is also optional.  If the spin is unknown (or zero), spinTuple can be set to None.
        #"""
        #spaceGroup = SpaceGroups.GetSpaceGroup(spaceGroupInt)
        
        #unitcell = Cell(spaceGroup, 0,0,0, a, b, c, alpha, gamma, beta)
        
        
        #for i in range(len(atomData)):
            ##Since I am adding anisotropy later, some things(files) will not have it
            ##anisotropy = None
            ##try:
            ##    anisotropy = (atomData[i][5], atomData[i][6], atomData[i][7])
            ##except IndexError:
            ##    print "Anisotropy not included!"
            ##if anisotropy == None:
            ##    anisotropy = (0,0,0)#Can't be None
            #anisotropy = (atomData[i][5], atomData[i][6], atomData[i][7])
            #if len(atomData[i]) > 11 and atomData[i][11] != None:#The entry can be left out or made None
                ##unitcell.generateAtoms((float(atomData[i][2]), float(atomData[i][3]), float(atomData[i][4])), atomData[i][0], atomData[i][1], atomData[i][9], anisotropy = anisotropy, spinMagnitude = atomData[i][8], elementObj = atomData[i][10])
                #unitcell.generateAtoms(elementObj = atomData[11], position = (atomData[2], atomData[3], atomData[4]),
                                       #massNumber = atomData[1], valence = atomData[9],
                                       #anisotropy = (atomData[5],atomData[6],atomData[7]),
                                       #spinMagnitude = atomData[8], description = "", spin = None, rgb = None, radius = None)
            #else:
                #elObj = getElObj(atomData[i][0])#Throws exception if symbol doesn't exist
                ##unitcell.generateAtoms((float(atomData[i][2]), float(atomData[i][3]), float(atomData[i][4])), atomData[i][0], atomData[i][1], atomData[i][9], anisotropy = anisotropy, spinMagnitude = atomData[i][8], elementObj = elObj)
                    
                    
        
        ##Create a Magnetic Cell
        ##self.MagCell = MagneticCell(unitcell, magNa, magNb, magNc, spaceGroup)
        ##Using Cutoff Cell as MagCell
        #self.MagCell = MagneticCell(unitcell, cutNa, cutNb, cutNc, spaceGroup)
        
        ##Regenerate Bonds as well
        #self.changeBonds(self.bondTable.data, notifyGUI)
    
    
    def openCif(self, filename):
        """Reads in a .cif file."""
        cf = CifFile.ReadCif(filename)
        
        #Assuming all data is in one outer block like NIST examples:
        data = cf[cf.keys()[0]]
        
        #Create a Crystollographic Unit Cell
        a = data['_cell_length_a']
        b = data['_cell_length_b']
        c = data['_cell_length_c']
        
        alpha = data['_cell_angle_alpha']
        gamma = data['_cell_angle_gamma']
        beta = data['_cell_angle_beta']
        
        spaceGroupInt = int(data['_symmetry_Int_Tables_number'])
        spaceGroup = SpaceGroups.GetSpaceGroup(spaceGroupInt)
        
        unitcell = Cell(spaceGroup, 0,0,0, a, b, c, alpha, gamma, beta)
        
        atomLabels = data['_atom_site_label']
        atomSymbol = data['_atom_site_type_symbol']
        xPositions = data['_atom_site_fract_x']
        yPositions = data['_atom_site_fract_y']
        zPositions = data['_atom_site_fract_z']
        
        atoms = [] #for the cell window
        for i in range(len(atomLabels)):
            #unitcell.generateAtoms((float(xPositions[i]), float(yPositions[i]), float(zPositions[i])), atomLabels[i])

            aData = [atomLabels[i], 0, float(xPositions[i]), float(yPositions[i]), float(zPositions[i])]
            #--Added to atomData: single ion anisotropy, spin magnitude, valence
            aData.append(0.0)#Dx
            aData.append(0.0)#Dy
            aData.append(0.0)#Dz
            aData.append(1)#Spin Magnitude
            aData.append('')#valence
            #-------------------------------------------------------------------
            atoms.append(aData)
            
            self.atomTable.SetValue(i, 0, atomLabels[i])
            self.atomTable.SetValue(i, 2, xPositions[i])
            self.atomTable.SetValue(i, 3, yPositions[i])
            self.atomTable.SetValue(i, 4, zPositions[i])
    
        #Create a Magnetic Cell
        self.MagCell = MagneticCell(unitcell, 1,1,1, spaceGroup)


        Na = 1  #Cif files only contain 1 unit cell
        Nb = 1
        Nc = 1
        
        #self.cellChange(spaceGroupInt, a, b, c, alpha, beta, gamma, magNa = Na, magNb = Nb, magNc = Nc, cutNa = Na, cutNb = Nb, cutNc = Nc, atomData = atoms)
        self.updateCell(spaceGroupInt, a, b, c, alpha, beta, gamma, magNa = Na, magNb = Nb, magNc = Nc, cutNa = Na, cutNb = Nb, cutNc = Nc, atomData = atoms)
        self.refreshGUI()
        
        
        #send signal to the cell window to show the info that has been loaded and to vtkWindow to draw it
        n =  self.atomTable.GetNumberRows()
        for i in range(n):
            print self.atomTable.GetValue(i, 0)
        send(signal = "File Load", sender = "Session", spaceGroup = spaceGroupInt, a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma, magNa = Na, magNb = Nb, magNc = Nc, cutNa = Na, cutNb = Nb, cutNc = Nc)
        
        
    def openFullProfFile(self, filename):
        """This will open a FullProf Studio file (.fst)"""
        handle = open(filename)
        lines = handle.readlines()
        handle.close()
        atoms = []
        bonds = []
        conns = []
        for line in lines:
            if line[0:4] == "CELL":
                #format of line: CELL a b c alpha beta gamma
                vals = line.split()
                print vals
                a = float(vals[1])
                b = float(vals[2])
                c = float(vals[3])
                alpha = float(vals[4])
                gamma = float(vals[5])
                beta = float(vals[6])
            elif line[0:6] == "SPACEG":
                #this is the space group in Hermann-Mauguin notation.
                hm_spacegroup = (line[6:]).strip().upper()
                space_group = GetSpaceGroup(hm_spacegroup)
            elif line[0:3] == "BOX":
                #Format: xmin xmax ymin ymax zmin zmax
                #In this program, however, xmin, ymin, zmin = 0,0,0 always.
                vals = line.split()
                a_diff = float(vals[2]) - float(vals[1])
                b_diff = float(vals[4]) - float(vals[3])
                c_diff = float(vals[6]) - float(vals[5])
                a_cutoff = int(a_diff)
                b_cutoff = int(b_diff)
                c_cutoff = int(c_diff)
                if a_diff - a_cutoff > 0:
                    a_cutoff += 1
                if b_diff - b_cutoff > 0:
                    b_cutoff += 1
                if c_diff - c_cutoff > 0:
                    c_cutoff += 1
            elif line[0:4] == "ATOM":
                vals = line.split()
                label = vals[1]
                symbol = vals[2]
                a_coord = float(vals[3])
                b_coord = float(vals[4])
                c_coord = float(vals[5])
                position = (a_coord, b_coord, c_coord)
                #Get the radius which is right after the word "RADIUS"
                for i in range(len(vals)):
                    if vals[i] == "RADIUS":
                        radius = float(vals[i+1])
                        break
                else:
                    radius = None
                #Get the color which is right after the word "COLOR"
                for i in range(len(vals)):
                    if vals[i] == "COLOR":
                        color = [float(vals[i+1]), float(vals[i+2]), float(vals[i+3])]
                        break
                else:
                    color = None
                #atomData format (each line):
                #label massNum aPos bPos cPos anisotropy_a anisotropy_b anistropy_c spin valence
                atoms.append([label, symbol, position, radius, color])
            elif line[0:4] == "BOND":
                #Format: BOND label1 label2 min_dist max_dist RADIUS rad COLOR r g b t
                #The color and radius need not be there and will be ignored for now since
                #the color and radius of bonds is hardcoded in right now.
                vals = line.split()
                bonds.append([vals[1], vals[2], vals[3], vals[4]])
            elif line[0:4] == "CONN":
                #Format: BOND symbol1 symbol2 min_dist max_dist RADIUS rad COLOR r g b t
                #The color and radius need not be there and will be ignored for now since
                #the color and radius of bonds is hardcoded in right now.
                vals = line.split()
                conns.append([vals[1], vals[2], vals[3], vals[4]])
                
        
        self.newCell(space_group.number, a, b, c, alpha, beta, gamma, 1, 1, 1,
                     a_cutoff, b_cutoff, c_cutoff)
        
        for atom in atoms:
            #FPStudio does not seem to support isotopes
            massNum = None
            self.addAtom(atom[1], atom[2], massNum = massNum, radius = atom[3], rgb = atom[4])
            
        for bond in bonds:
            self.createBonds(label1 = bonds[0], label2 = bonds[1],
                             minDist = bonds[2], maxDist = bonds[3])
        for conn in conns:
            self.createBonds(symbol1 = conns[0], symbol2 = conns[1],
                             minDist = conns[2], maxDist = conns[3])
            
        self.refreshGUI()
        #self.cellChange(space_group.number, a, b, c, alpha, beta, gamma, magNa = 1, magNb = 1, magNc = 1, cutNa = a_cutoff, cutNb = b_cutoff, cutNc = c_cutoff, atomData = atoms)
        #self.updateCell(space_group.number, a, b, c, alpha, beta, gamma, magNa = 1, magNb = 1, magNc = 1, cutNa = a_cutoff, cutNb = b_cutoff, cutNc = c_cutoff, atomData = atoms)
        #self.refreshGUI()
        
        #send signal to the cell window to show the info that has been loaded and to vtkWindow to draw it
        send(signal = "File Load", sender = "Session",
             spaceGroup = space_group.number, a = a, b = b, c = c,
             alpha = alpha, beta = beta, gamma = gamma, magNa = a_cutoff,
             magNb = b_cutoff, magNc = c_cutoff, cutNa = a_cutoff,
             cutNb = b_cutoff, cutNc = c_cutoff)
        
               
            #TODO:  use these values extracted.  You could combine the three file opening functions.
            #Each function would have to extract values form it's format and then a single function
            #could be used for all three to construct the model from the extracted values.e
                

    def createBonds(self, minDist, maxDist, label1 = None, label2 = None,
                    symbol1 = None, symbol2 = None):
        """Creates bonds between atoms at the right distance apart.
        
        Either symbol1 and symbol 2 must be supplied, or label1 and label2.  A
        bond will be created between all atoms whose seperation is between
        minDist and maxDist if one has label1 and the other label2, or one has
        element symbol1, and the other symbol2."""
        #This is not as efficient as it could be, since symmetry equivalent
        #bonds will be checked.  Only atoms in the first cell and surrounding
        #cell need to be checked.
        if symbol1 != None and symbol2 != None:
            for a1 in self.MagCell.getAllAtoms():
                for a2 in self.MagCell.getAllAtoms():
                    if a1.getSymbol() == symbol1 and a2.getSymbol() == symbol2:
                        self.addBond(atom1 = a1, atom2 = a2)
                        return
        if label1 != None and label2 != None:
            for a1 in self.MagCell.getAllAtoms():
                for a2 in self.MagCell.getAllAtoms():
                    if a1.getDescription() == symbol1:
                        if a2.getDescription() == symbol2:
                            self.addBond(atom1 = a1, atom2 = a2)
                            return
                        
                
    
    def saveSessionToXML(self, filename):
        """Saves all the information needed to reconstruct the model in an xml file."""
        xmlStr = self.createXMLStr()
        
        #Write to the file
        #xml.dom.ext.PrettyPrint(doc, open(filename, 'w'))
        xmlFile = open(filename, 'w')
        xmlFile.write(xmlStr)
        xmlFile.close()
        
    def createXMLStr(self):
        #Create the document
        doc = xml.dom.minidom.Document()
        
        #Create the main element
        main_element = doc.createElement('Spinwaves_Session')
        doc.appendChild(main_element)
        
        #Currently this will save cell info from the last time the 
        #generate button was pressed, but atom and bond info that is
        #currently displayed in the tables
        
        #Add Unit Cell Attributes
        unitCell = doc.createElement('unit_cell')
        main_element.appendChild(unitCell)
        unitCell.setAttribute('space_group', str(self.MagCell.getAllUnitCells()[0].getSpaceGroup().number))
        unitCell.setAttribute('a', str(self.MagCell.getAllUnitCells()[0].getA()))
        unitCell.setAttribute('b', str(self.MagCell.getAllUnitCells()[0].getB()))
        unitCell.setAttribute('c', str(self.MagCell.getAllUnitCells()[0].getC()))
        unitCell.setAttribute('alpha', str(self.MagCell.getAllUnitCells()[0].getAlpha()))
        unitCell.setAttribute('beta', str(self.MagCell.getAllUnitCells()[0].getBeta()))
        unitCell.setAttribute('gamma', str(self.MagCell.getAllUnitCells()[0].getGamma()))
        
        #Add Magnetic Unit Cell Attributes
        magCell = doc.createElement('magnetic_cell')
        main_element.appendChild(magCell)
        magCell.setAttribute('Na', str(self.MagCell.getNa()))
        magCell.setAttribute('Nb', str(self.MagCell.getNb()))
        magCell.setAttribute('Nc', str(self.MagCell.getNc()))
        
        #Cutoff Cell Not yet being used
        
        
        #Add Atoms
        atomsElement = doc.createElement('atoms')
        main_element.appendChild(atomsElement)
        for i in range(self.atomTable.GetNumberRows()):
            atomElement = doc.createElement('atom')
            atomElement.setAttribute('Name', str(self.atomTable.GetValue(i, 0)))
            atomElement.setAttribute('Atomic_Number', str(self.atomTable.GetValue(i, 1)))
            atomElement.setAttribute('valence', self.atomTable.GetValue(i, 2))
            atomElement.setAttribute('x', str(self.atomTable.GetValue(i, 3)))
            atomElement.setAttribute('y', str(self.atomTable.GetValue(i, 4)))
            atomElement.setAttribute('z', str(self.atomTable.GetValue(i, 5)))
            atomElement.setAttribute('SIA_Dx', str(self.atomTable.GetValue(i, 6)))
            atomElement.setAttribute('SIA_Dy', str(self.atomTable.GetValue(i, 7)))
            atomElement.setAttribute('SIA_Dz', str(self.atomTable.GetValue(i, 8)))
            atomElement.setAttribute('spin_magnitude', str(self.atomTable.GetValue(i, 9)))
            atomsElement.appendChild(atomElement)
            
          
        #Add Bonds
        bondsElement = doc.createElement('bonds')
        main_element.appendChild(bondsElement)
        for i in range(self.bondTable.GetNumberRows()):
            bondElement = doc.createElement('bond')
            bondElement.setAttribute('Atom1_Number', str(self.bondTable.GetValue(i, 0)))
            bondElement.setAttribute('Atom1_Na', str(self.bondTable.GetValue(i, 1)))
            bondElement.setAttribute('Atom1_Nb', str(self.bondTable.GetValue(i, 2)))
            bondElement.setAttribute('Atom1_Nc', str(self.bondTable.GetValue(i, 3)))
            bondElement.setAttribute('Atom2_Number', str(self.bondTable.GetValue(i, 4)))
            bondElement.setAttribute('Atom2_Na', str(self.bondTable.GetValue(i, 5)))
            bondElement.setAttribute('Atom2_Nb', str(self.bondTable.GetValue(i, 6)))
            bondElement.setAttribute('Atom2_Nc', str(self.bondTable.GetValue(i, 7)))
            bondElement.setAttribute('On', str(self.bondTable.GetValue(i, 9)))
            bondsElement.appendChild(bondElement)
          
            
            #Right now I allow bond creation without Jmatrix         
            matrix = self.bondTable.GetActualValue(i, 8)
            print matrix
            #if matrix != '':
            jMatrix = doc.createElement('jMatrix')
            for i in range(3):
                for j in range(3):
                    parameter = doc.createElement('Parameter')
                    parameter.setAttribute('position', 'j'+str(i+1)+str(j+1))
                    parameter.setAttribute('name', matrix[i][j].getName())
                    parameter.setAttribute('fit', str(matrix[i][j].fit))
                    parameter.setAttribute('value', str(matrix[i][j].value))
                    parameter.setAttribute('min', matrix[i][j].min)
                    parameter.setAttribute('max', matrix[i][j].max)
                    parameter.setAttribute('default', str(matrix[i][j].default))
                    parameter.setAttribute('tiedTo', str(matrix[i][j].tied))
                    jMatrix.appendChild(parameter)
            bondElement.appendChild(jMatrix)
            
        #Write to screen
        #xml.dom.ext.PrettyPrint(doc)
        xmlStr = doc.toprettyxml("    ")
        return xmlStr
       
    #Not used by me, but may be useful for other people, should add this to GUI  
    def export(self, filename):
        """Exports the bond information to a file."""
        size = 10
        
        initialTime = time.clock()
        
        file = open(filename, 'w')
#        file.write("#Atoms\n#Number X Y Z\n")
#        for i in range(size):
#            for j in range(size):
#                for k in range(size):
#                    for atom in self.getCutoffCell().getAllAtoms():
#                        num = atom.getIndexNumber()
#                        pos = atom.getPosition()
#                        x = pos[0] + (i * self.getCutoffCell().getNa())
#                        y = pos[1] + (j * self.getCutoffCell().getNa())
#                        z = pos[2] + (k * self.getCutoffCell().getNa())
#                        line = str(num) + " " + str(x) + " " + str(y) + " " + str(z)
#                        file.write(line + "\n")
        #Right now this will be extremely slow (8 nested for loops)       
        class SimpleBond():
            def __init__(self, pos1, pos2, jMatrix):
                self.pos1 = pos1
                self.pos2 = pos2
                self.jMatrix = jMatrix
                
            def sameBond(self, bond2):
                if self.pos1 == bond2.pos1 or self.pos1 == bond2.pos2:
                    if self.pos2 == bond2.pos2 or self.pos2 == bond2.pos1:
                        return True
                return False
        
        class SimpleBondList():
            def __init__(self):
                self.list = []
            
            def addBond(self, bond):
                if not self.containsBond(bond):
                    self.list.append(bond)
                    
            def containsBond(self, bond):
                for eachBond in self.list:
                    if eachBond.sameBond(bond):
                        return True
                return False
        
        Na = self.getCutoffCell().getNa()
        Nb = self.getCutoffCell().getNb()
        Nc = self.getCutoffCell().getNc()
        
        
        def contains(list, element):
            for item in list:
                if (item == element).all():
                    return True
            return False
        
        def indexOf(list, item):
            for i in range(len(list)):
                if item.all() == list[i].all():
                    return i
            return -1
        
        
        matrices = []
        for bond in self.getCutoffCell().getBonds():
#           pos1 = bond.getAtom1().getPosition()
#           pos2 = bond.getAtom2().getPosition()
            jMat = bond.getJMatrix()
#            count = matrices.count(jMat)
            if not contains(matrices, jMat):
                matrices.append(jMat)
        
        
        simpleCellBonds = []
        for bond in self.getCutoffCell().getBonds():
            pos1 = bond.getAtom1().getPosition()
            pos2 = bond.getAtom2().getPosition()
            jMat = bond.getJMatrix()
            newBond = SimpleBond(pos1, pos2, indexOf(matrices,jMat))
            simpleCellBonds.append(newBond)
        
        
        simpleBonds = SimpleBondList()
        for bond in simpleCellBonds:
            pos1 = bond.pos1
            pos2 = bond.pos2
            jMatInt = bond.jMatrix
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for a in range(Na):
                            for b in range(Nb):
                                for c in range(Nc):
                                    x1 = pos1[0] + a + (Na * i)
                                    y1 = pos1[1] + b + (Nb * j)
                                    z1 = pos1[2] + c + (Nc * k)
                                    
                                    x2 = pos2[0] + a + (Na * i)
                                    y2 = pos2[1] + b + (Nb * j)
                                    z2 = pos2[2] + c + (Nc * k)    
                                    bond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMatInt )
                                    simpleBonds.addBond(bond)                              

        #Pick out bonds that link first cutoff cell and another
        interCutoffBonds = []
        for eachBond in simpleBonds.list:
            #Check if one bond is in the first cutoff cell
            if eachBond.pos1[0] < Na or eachBond.pos2[0] < Na: #x
                if eachBond.pos1[1] < Nb or eachBond.pos2[1] < Nb: #y
                    if eachBond.pos1[2] < Nc or eachBond.pos2[2] < Nc: #z
                        #check if the second bond is not in the first cutoff cell
                        if (not eachBond.pos1[0] < Na) or (not eachBond.pos2[0] < Na): #x
                            if (not eachBond.pos1[1] < Nb) or (not eachBond.pos2[1] < Nb): #y
                                if (not eachBond.pos1[2] < Nc) or (not eachBond.pos2[2] < Nc): #z
                                    interCutoffBonds.append(eachBond)
                                    
        
        
        
#        finalBondList = []
        #Translate all bonds within the cutoff cell

            
        file.write("#J Matrices\n#Number J11 J12 J13 J21 J22 J23 J31 J32 J33\n")
        for i in range(len(matrices)):
            jMat = matrices[i]
            jStr = str(i) + " " + str(jMat[0][0]) + " " + str(jMat[0][1]) + " " + str(jMat[0][2]) + " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            file.write(jStr + "\n")
        
        
        file.write("#Bonds\n#X1 Y1 Z1 X2 Y2 Z2 J\n")
        for bond in simpleCellBonds:   
            pos1 = bond.pos1
            pos2 = bond.pos2
#            jStr += " " 
#            jStr += str(jMat[0][1]) 
#            jStr += " " + str(jMat[0][2]) 
#            jStr += " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            for i in range(size):
                for j in range(size):
                    for k in range(size):
                        x1 = pos1[0] + (Na * i)
                        y1 = pos1[1] + (Nb * j)
                        z1 = pos1[2] + (Nc * k)
                        
                        x2 = pos2[0] + (Na * i)
                        y2 = pos2[1] + (Nb * j)
                        z2 = pos2[2] + (Nc * k)    
#                        smplBond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMat )
 #                       finalBondList.append(smplBond)
                        pos1Str = str(x1) + " " + str(y1) + " " + str(z1)
                        pos2Str = str(x2) + " " + str(y2) + " " + str(z2)
                        jStr = str(bond.jMatrix)
#There should always be a jMatrix                        if jMat != None:
                        file.write(pos1Str + " " + pos2Str + " " + jStr + "\n")
#                    else:
#                        file.write(pos1Str + " " + pos2Str + "\n")
        
        for smplBond in interCutoffBonds:
            pos1 = smplBond.pos1
            pos2 = smplBond.pos2
            jMat = smplBond.jMatrix
            jStr = str(jMat)
#            jStr = str(jMat[0][0]) + " " + str(jMat[0][1]) + " " + str(jMat[0][2]) + " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            aDisp = abs(pos1[0] - pos2[0])
            bDisp = abs(pos1[1] - pos2[1])
            cDisp = abs(pos1[2] - pos2[2])
#            if pos1[0] > pos2[0]:
#                aDisp = pos1[0]
#            else:
#                aDisp = pos2[0]
#            if pos1[1] > pos2[1]:
#                bDisp = pos1[1]
#            else:
#                bDisp = pos2[1]
#            if pos1[2] > pos2[2]:
#                cDisp = pos1[2]
#            else:
#                cDisp = pos2[2]
            for i in range(size - aDisp):
                for j in range(size - bDisp):
                    for k in range(size - cDisp):
                        x1 = pos1[0] + (Na * i)
                        y1 = pos1[1] + (Nb * j)
                        z1 = pos1[2] + (Nc * k)
                        
                        x2 = pos2[0] + (Na * i)
                        y2 = pos2[1] + (Nb * j)
                        z2 = pos2[2] + (Nc * k)    
#                        smplBond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMat )
 #                       finalBondList.append(smplBond)
                        pos1Str = str(x1) + " " + str(y1) + " " + str(z1)
                        pos2Str = str(x2) + " " + str(y2) + " " + str(z2)
#There should always be a jMatrix                        if jMat != None:
                        file.write(pos1Str + " " + pos2Str + " " + jStr + "\n")
        
        #Check for reapeats in finalBond List just for testing
#        def isRepeat(finalBondList):
#            for i in range(0, len(finalBondList)):
#                for j in range(i + 1, len(finalBondList)):
#                    if finalBondList[i].sameBond(finalBondList[j]):
#                        return True
#            return False
#        
#        if isRepeat(finalBondList):
#            print "There is a repeat!"
#        else:
#            print "NO repeats!"
        
        

        
        file.close()
        seconds = time.clock() - initialTime
        minutes = int(seconds)/60
        seconds -= (minutes*60)
        hours = minutes/60
        minutes -= hours*60
        print "Done\nTime:", hours, "hours", minutes, "minutes", seconds, "seconds" 
        
       
        
    #This one orgnizes it by atom, not bond
    
    
    def exportForMonteCarlo(self, filename, size):
        """Exports the bond information to a file in a format more useful for
        the Simulated Annealing.

        #AtomNumber AtomPosition(X Y Z) Anisotropy(X Y Z) OtherIndex Jmatrix OtherIndex Jmatrix..."""
        #DO NOT MAKE THE SIZE LESS THAN 3.  The implementation depends on it.  Particulary the
        #validBond method, but the dispersion and cross section code is also built around this.
        if size < 3:
            size = 3
            
        timer = Timer()
        matrices, allAtoms = self.Export_Aux(size)  
        timer.printTime()
        print "number of atoms: ", len(allAtoms), "\n writing to disk..."
        
        file = open(filename, 'w')
        
        Na = self.getCutoffCell().getNa()
        Nb = self.getCutoffCell().getNb()
        Nc = self.getCutoffCell().getNc()
            
        #Can add flag in here if the coordinates are <= Na, Nb, Nc
        #(if it's in the cutoff cell) for the spinwave calculation
        def inInteractionCellStr(atoms, atom):
            """Used for output to create an "X" if the atom is in the first interaction
            Cell or "O" if not.  This is the actual smallest interaction cell, not the
            cutoff cell created by the user.  An atom is in the first cutoff cell if it
            is either in the first crystallographic unit cell or if it bonds with an atom
            that is.  The 'first crystallographic unit cell' will not be the cell
            at (0,0,0), but rather the corresponding cell in the cutoff cell at
            (1,1,1) (measured in cutoff cells, not unit cells).  The desired
            crystallographic cell is therefore at (Na, Nb, Nc).  This is to ensure
            that the cell is completely surrounded and therefore no interactions
            will be left out.  If the cell at (0,0,0) were used, interactions
            in the negative direction would not be included."""
            
            def inDesiredCell(atom):
                if atom.pos[0] >= Na and atom.pos[0] < (Na + 1):
                    if atom.pos[1] >= Nb and atom.pos[1] < (Nb + 1):
                        if atom.pos[2] >= Nc and atom.pos[2] < (Nc + 1):
                            return True
                return False
            
            #First check if the atom is in the first crystallographic cell
            if inDesiredCell(atom):
                return "X"
            #If not, check if it bonds to an atom that is
            for i in range(len(atom.interactions)):
                if inDesiredCell(atoms[atom.interactions[i][0]]):
                    return "X"
            
            for interaction in atom.interCellInteractions:
                interactingAtom = atoms[interaction[0]]
                if inDesiredCell(interactingAtom):
                    return "X"
            
            return "O"
        
        
        #print the size of the cutoff(interaction) cell
        #This is necessary now that the first unit cell is not at (0,0,0), but
        #at (Na, Nb, Nc)
        file.write("#Interaction Cell Dimensions\n#Na Nb Nc\n")
        file.write(str(Na) + " " + str(Nb) + " " + str(Nc) + "\n")
        
        #Write the matrix list to the file
        file.write("#J Matrices\n#Number J11 J12 J13 J21 J22 J23 J31 J32 J33\n")
        for i in range(len(matrices)):
            jMat = matrices[i]
            #Old method for J matrices of floats
            #jStr = str(i) + " " + str(jMat[0][0]) + " " + str(jMat[0][1]) + " " + str(jMat[0][2]) + " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            jStr = str(i) + " " + str(jMat[0][0].default)
            jStr += " " + str(jMat[0][1].default) 
            jStr += " " + str(jMat[0][2].default)
            jStr += " " + str(jMat[1][0].default) 
            jStr += " " + str(jMat[1][1].default) 
            jStr += " " + str(jMat[1][2].default)
            jStr += " " + str(jMat[2][0].default) 
            jStr += " " + str(jMat[2][1].default) 
            jStr += " " + str(jMat[2][2].default)
            file.write(jStr + "\n")
        
        
        
        #print out the simple atom list
        file.write("#AtomNumber Symbol CellNumber MassNumber  Valence InFirstInteractionCell AtomPosition(X Y Z) Anisotropy(X Y Z) SpinMagnitude OtherIndex Jmatrix OtherIndex Jmatrix...\n")
        for atomIndex in range(len(allAtoms)):
            atom = allAtoms[atomIndex]
            atomStr = str(atomIndex) + " " + str(atom.label) + " " + str(atom.cellNum)
            atomStr += " " + str(atom.massNum) + " " + str(atom.valence)
            if atom.valence == "":
                atomStr += "_"
            atomStr += " " + inInteractionCellStr(allAtoms, atom)
            atomStr += " " + str(atom.pos[0]) + " " + str(atom.pos[1]) + " " + str(atom.pos[2])
            atomStr += " " + str(atom.anisotropy[0]) + " " + str(atom.anisotropy[1]) + " " + str(atom.anisotropy[2])
            atomStr += " " + str(atom.spinMag)
            for interaction in atom.interactions:
                otherAtom = interaction[0]
                jMat = interaction[1]
                atomStr += " " + str(otherAtom)
                atomStr += " " + str(jMat)
            for interaction in atom.interCellInteractions:#And again for inter-cell
                otherAtom = interaction[0]
                jMat = interaction[1]
                atomStr += " " + str(otherAtom)
                atomStr += " " + str(jMat)
            file.write(atomStr + "\n")
        
        file.close()
        print "done"
        timer.printTime()
            
      
    def Export_Aux(self, size):
        """The section of export that does not write to the file, but organizes
        the atoms in a way that is useful for the monte carlo simulation has
        been moved to this function, so that the simulation can be run without
        files."""
        class SimpleBond():
            def __init__(self, pos1, pos2, jMatrix, anisotropy1, anisotropy2, spinMag1, spinMag2):
                self.pos1 = pos1
                self.anisotropy1 = anisotropy1
                self.pos2 = pos2
                self.anisotropy2 = anisotropy2
                self.jMatrix = jMatrix
                self.spinMag1 = spinMag1
                self.spinMag2 = spinMag2
                
            def sameBond(self, bond2):
                if self.pos1 == bond2.pos1 or self.pos1 == bond2.pos2:
                    if self.pos2 == bond2.pos2 or self.pos2 == bond2.pos1:
                        return True
                return False
        
        Na = self.getCutoffCell().getNa()
        Nb = self.getCutoffCell().getNb()
        Nc = self.getCutoffCell().getNc()
        
        class SimpleAtom():
            def __init__(self, pos, anisotropy, spinMag, label, massNum, cellNum, valence):
                self.anisotropy = anisotropy
                self.pos = pos
#                self.cellPos = []
#                self.cellPosX = int(pos[0])/Na
#                self.cellPosY = int(pos[1])/Nb
#                self.cellPosZ = int(pos[2])/Nc
                self.interactions = []
                self.interCellInteractions = []#these must be translated with different logic and therefore must be kept separate
                #self.interactions[position of other atom] = j number
                self.spinMag = spinMag
                self.label = label
                self.massNum = massNum
                self.cellNum = cellNum
                #valence is an int, not str
                #self.valence = string.join(valence.split(), "")#just in case, get rid of whitespace
                self.valence = valence
                
            #might want to change this to position later when all atoms wont be created in same list
            def addInteraction(self, atom2, jMat):
#                self.interactions[atom2] = jMat
                self.interactions.append((atom2, jMat))
                
            def addInterCellInteraction(self, atom2, jMat, direction):
                """Direction is in form (bool, bool, bool) for (x,y,z)"""
                #Check for repeats (necessary for method used to translate these bonds)
                for interaction in self.interCellInteractions:
                    if interaction[0] == atom2 and interaction[1] == jMat:
                        return #contains this interaction already
                self.interCellInteractions.append((atom2, jMat, direction))
            
            def __eq__(self, other):
                if self.pos[0] == other.pos[0]:
                    if self.pos[1] == other.pos[1]:
                        if self.pos[2] == other.pos[2]:
                            return True
                return False
        
        
        class SimpleBondList():
            def __init__(self):
                self.list = []
                
            def addBond(self, bond):
                if not self.containsBond(bond):
                    self.list.append(bond)
#                else:
#                    print "Duplicate Bonds!" #should not get here
                    
            def containsBond(self, bond):
                for eachBond in self.list:
                    if eachBond.sameBond(bond):
                        return True
                return False
        
        
        
        def contains(list, element):
            for item in list:
                if (item == element).all():
                    return True
            return False
        
        def atomListContains(list, element):
            for item in list:
                if item == element:
                    return True
            return False
        
        def indexOf(list, item):
            for i in range(len(list)):
                if (item == list[i]).all():
                    return i
            return -1
        
        def translateToFirstCutoffCell(pos):
            """Translates a position back to the first Cutoff cell."""
            x = pos[0]
            y = pos[1]
            z = pos[2]
            
            while x >= Na:
                x = x - Na
                
            while y >= Nb:
                y = y - Nb
                
            while z >= Nc:
                z = z - Nc
                
            return (x,y,z)


        
        #Create list of matrices
        matrices = []
        for bond in self.getCutoffCell().getBonds():
#           pos1 = bond.getAtom1().getPosition()
#           pos2 = bond.getAtom2().getPosition()
            jMat = bond.getJMatrix()
#            count = matrices.count(jMat)
            if not contains(matrices, jMat):
                matrices.append(jMat)
        
        
        #create simple bonds within cutoff cell
        simpleCellBonds = []
        for bond in self.getCutoffCell().getBonds():
            pos1 = bond.getAtom1().getPosition()
            anisotropy1 = bond.getAtom1().getAnisotropy()
            spin1 = bond.getAtom1().getSpinMagnitude()
#            print anisotropy1
            pos2 = bond.getAtom2().getPosition()
            anisotropy2 = bond.getAtom2().getAnisotropy()
            spin2 = bond.getAtom2().getSpinMagnitude()
#            print anisotropy2
#            time.sleep(.1)
            jMat = bond.getJMatrix()
            newBond = SimpleBond(pos1, pos2, indexOf(matrices,jMat), anisotropy1, anisotropy2, spin1, spin2)
            simpleCellBonds.append(newBond)
        
        
        def PosInFirstCutoff(pos):
            return (pos[0] < Na and pos[1] < Nb and pos[2] < Nc)
        
        
        cellAtoms = []
        cellBonds = SimpleBondList()
        interCellBonds = SimpleBondList()#bonds between cells cannot be translated
        #with the same index arithmetic because of edges
        for bond in simpleCellBonds:
            pos1 = bond.pos1
            anisotropy1 = bond.anisotropy1
            pos2 = bond.pos2
            anisotropy2 = bond.anisotropy2
            spin1 = bond.spinMag1
            spin2 = bond.spinMag2
            #--Adding labels, cell number, valence, and mass number----
            a1 = self.MagCell.atomAtPosition(pos1)
            a2 = self.MagCell.atomAtPosition(pos2)
            #label1 = a1.description
            #label2 = a2.description
            label1 = a1.getSymbol()
            label2 = a2.getSymbol()
            #cellNum1 = a1.getIndexNumber()
            #cellNum2 = a2.getIndexNumber()
            cellNum1 = a1.getIDNum()
            cellNum2 = a2.getIDNum()
            massNum1 = a1.getMassNum()
            massNum2 = a2.getMassNum()
            valence1 = a1.getValence()
            valence2 = a2.getValence()
            #---------------------------------------------------
            jMatInt = bond.jMatrix
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for a in range(Na):
                            for b in range(Nb):
                                for c in range(Nc):
                                    x1 = pos1[0] + a + (Na * i)
                                    y1 = pos1[1] + b + (Nb * j)
                                    z1 = pos1[2] + c + (Nc * k)
                                    
                                    x2 = pos2[0] + a + (Na * i)
                                    y2 = pos2[1] + b + (Nb * j)
                                    z2 = pos2[2] + c + (Nc * k)  
                                    newPos1 = (x1,y1,z1)
                                    newPos2 = (x2,y2,z2)
                                    
                                    #It is possible for an atom to be bonded to an atom which was
                                    #translated from a non-bonded atom in the original unit cell.
                                    #This non-bonded atom which translates into bonded atom(s) must
                                    #be included in the cellAtoms list, or it will not be translated
                                    #through the allAtoms list to create the other atoms which are
                                    #bonded, and hence, must be included.
                                    
                                    if PosInFirstCutoff(newPos1) or PosInFirstCutoff(newPos2):  
                                        if PosInFirstCutoff(newPos1) and PosInFirstCutoff(newPos2):
                                            #Both are in first cutoff
                                            #Add the atoms to the list of atoms within the first cell
                                            newAtom1 = SimpleAtom(newPos1, anisotropy1, spin1, label1, massNum1, cellNum1, valence1)
                                            if not atomListContains(cellAtoms, newAtom1):
                                                cellAtoms.append(newAtom1)
                                            newAtom2 = SimpleAtom(newPos2, anisotropy2, spin2, label2, massNum2, cellNum2, valence2)
                                            if not atomListContains(cellAtoms, newAtom2):
                                                cellAtoms.append(newAtom2)
                                            #Add the bond to bonds within the cell
                                            bond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMatInt , None, None, None, None)
                                            cellBonds.addBond(bond)
                                        else:#It is an inter-cellular bond
                                            bond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMatInt, None, None, None, None)
                                            interCellBonds.addBond(bond)
                                            #If the atom is in the first cutoff cell then it must be added and
                                            #translating the position will do nothing.  If it is not in the first cutoff
                                            #cell, then the corresponding atom in the first cutoff cell must be added
                                            #to create this one through translation
                                            transPos1 = translateToFirstCutoffCell(newPos1)
                                            transPos2 = translateToFirstCutoffCell(newPos2)
                                            
                                            newAtom1 = SimpleAtom(transPos1, anisotropy1, spin1, label1, massNum1, cellNum1, valence1)
                                            newAtom2 = SimpleAtom(transPos2, anisotropy2, spin2, label2, massNum2, cellNum2, valence2)
                                            if not atomListContains(cellAtoms, newAtom1):
                                                cellAtoms.append(newAtom1)
                                            if not atomListContains(cellAtoms, newAtom2):
                                                cellAtoms.append(newAtom2)
                                    
        print "after iterating through bonds..."
        print "atoms:"
        for atom in cellAtoms:
            print atom.pos
        print "cellBonds: "
        for b in cellBonds.list:
            print b.pos1, " , ", b.pos2
        print "interCelBbonds: "
        for b in interCellBonds.list:
            print b.pos1, " , ", b.pos2
        
        #symmetry equivalent bonds between unit cells will not be represented in
        #the cutoff cell if the cutoff cell is only one unit cell wide in any
        #dimension which would include these inter-cellular bonds
        if size > 1 and (Na == 1 or Nb == 1 or Nc == 1):
            for bond in simpleCellBonds:
                xyz = bond.pos1
                xyz2 = bond.pos2
                
                #one of the two atoms should be in the first unit cell
                if(xyz[0] < 1 and xyz[1] < 1 and xyz[2] < 1) or (xyz2[0] < 1 and xyz2[1] < 1 and xyz2[2] < 1):
                
                    for symop in self.MagCell.space_Group.iter_symops():
                    # operate on coordinates in non-shifted spacegroup
                        pos1 = symop(xyz)
                        pos2 = symop(xyz2)
                        
                        mask1 = numpy.logical_or(pos1 < 0.0, pos1 >= 1.0)
                        translation = numpy.floor(pos1[mask1])  #translates the first atom back to cell at (0,0,0)
                        pos1[mask1] -= translation
                        pos2[mask1] -= translation  #Uses same translation to translate other atom
                             
                        
                        #translate new Bond by 1 cell in each direction so all
                        #translations of intercellular bonds are represented.
                        
    
                        #iterate through each translation and check if there are atoms there that could
                        #be bonded; if so, add the bond
                        for i in range(0, 2): #translate in x direction (Na - Cell X position) times
                            for j in range(0, 2): #translate in y direction (Nb - Cell Y position) times
                                for k in range(0, 2): #translate in z direction (Nc - Cell Z position) times
                                    translatedPos1 = [i + pos1[0],j + pos1[1],k + pos1[2]]
                                    translatedPos2 = [i + pos2[0],j + pos2[1],k + pos2[2]]
                                    
                                    
                                    #Check if the bond crosses the border(s) of the dimension of only one unit cell.
                                    #Check if the bond exists in intercellular bonds, and if not, add it? redundant ^?
                                    #Then check if the aotm that is in the cutoff cell is represented in cellAtoms, and
                                    #if not, add it.
                                    
                                    #If the bond crosses a dimension of size 1 unit cell
                                    if ((Na == 1) and (int(translatedPos1[0]) != int(translatedPos2[0]) or translatedPos2[0] < 0)) or ((Nb == 1) and (int(translatedPos1[1]) != int(translatedPos2[1]) or translatedPos2[1] < 0)) or ((Nc == 1) and (int(translatedPos1[2]) != int(translatedPos2[2]) or translatedPos2[2] < 0)):
                                        
                                        #Add the atom in the cutoff Cell and add the bond to intercellular bonds
                                        print translatedPos1, translatedPos2
                                        atomObj1 = self.MagCell.atomAtPosition(translatedPos1)
                                        atomObj2 = self.MagCell.atomAtPosition(translatedPos2)
                                        if(atomObj1 != None):#Add the atom if it is in the cutoff cell
                                            newAtom1 = SimpleAtom(translatedPos1, atomObj1.anisotropy, atomObj1.spinMagnitude, atomObj1.description, atomObj1.getMassNum(), atomObj1.getIDNum(), atomObj1.valence)
                                        #Add atom if there is not already an atom at that position
                                        if not atomListContains(cellAtoms, newAtom1):
                                            cellAtoms.append(newAtom1)
                                        
                                        if(atomObj2 != None):#Add the atom if it is in the cutoff cell
                                            newAtom2 = SimpleAtom(translatedPos2, atomObj2.anisotropy, atomObj2.spinMagnitude, atomObj2.description, atomObj2.getMassNum(), atomObj2.getIDNum(), atomObj2.valence)
                                        #Add atom if there is not already an atom at that position
                                        if not atomListContains(cellAtoms, newAtom2):
                                            cellAtoms.append(newAtom2)
                                            
                                        #If one of the atoms are in the cutoff cell and both have positive coordinates, add the bond
                                        if(atomObj1 != None) or atomObj2 != None:
                                            if(translatedPos1[0] >= 0 and translatedPos1[1] >= 0 and translatedPos1[2] >= 0 and translatedPos2[0] >= 0 and translatedPos2[1] >= 0 and translatedPos2[2] >= 0):
                                                interCellBonds.addBond(SimpleBond(translatedPos1, translatedPos2, bond.jMatrix, None, None, None, None))
                                        
                                        
        print "after iterating through symops..."
        print "atoms:"
        for atom in cellAtoms:
            print atom.pos
        print "cellBonds: "
        for b in cellBonds.list:
            print b.pos1, " , ", b.pos2
        print "interCelBbonds: "
        for b in interCellBonds.list:
            print b.pos1, " , ", b.pos2
        allAtoms = []
        numAtomsPerCell = len(cellAtoms)
        
        print "atoms in cell: ", numAtomsPerCell
        
        for i in range(size):
            for j in range(size):
                for k in range(size):
                    for index in range(len(cellAtoms)):
                        pos = cellAtoms[index].pos
                        anisotropy = cellAtoms[index].anisotropy
                        x = pos[0] + (Na * i)
                        y = pos[1] + (Nb * j)
                        z = pos[2] + (Nc * k)
                        newAtom = SimpleAtom((x,y,z), anisotropy, cellAtoms[index].spinMag, cellAtoms[index].label, cellAtoms[index].massNum, cellAtoms[index].cellNum, cellAtoms[index].valence)
#                        print (len(allAtoms)%numAtomsPerCell == index)#just a check, should always be true
                        allAtoms.append(newAtom)

 
        #for atom in allAtoms:
        #    print atom.pos
        
        #Add bonds cellBonds to allAtoms (currently not most efficient way, but its a short list)
        for bond in cellBonds.list:
            #Make sure each position is not yet represented the easy but inefficient way
            #It's not actually that inefficient since it should only go through a fraction of the list
            pos1 = bond.pos1
            pos2 = bond.pos2
            pos1Index = -1
            pos2Index = -1
            for i in range(len(allAtoms)):
                currentPos = allAtoms[i].pos
                #if currentPos == pos1:
                if currentPos[0] == pos1[0] and currentPos[1] == pos1[1] and currentPos[2] == pos1[2]:
                    pos1Index = i
                    break
                
            for i in range(len(allAtoms)):
                currentPos = allAtoms[i].pos  
                #if currentPos == pos2:
                if currentPos[0] == pos2[0] and currentPos[1] == pos2[1] and currentPos[2] == pos2[2]:
                    pos2Index = i
                    break
            
            if pos1Index < 0 or pos2Index < 0:
                print "Atom list does not contain all atoms!"
                if pos1Index < 0:
                    print pos1, " missing"
                if pos2Index < 0:
                    print pos2, " missing"
                raise Exception("Export Failed")
            else:
                allAtoms[pos1Index].addInteraction(pos2Index, bond.jMatrix)
                allAtoms[pos2Index].addInteraction(pos1Index, bond.jMatrix)


        def bondDirection(pos1, pos2):
            xCell1 = int(pos1[0]/Na)  #Find the cutoff cell
            xCell2 = int(pos2[0]/Na)
            yCell1 = int(pos1[1]/Nb)
            yCell2 = int(pos2[1]/Nb)
            zCell1 = int(pos1[2]/Nc)
            zCell2 = int(pos2[2]/Nc)
            xShiftBool = (xCell1 != xCell2)
            yShiftBool = (yCell1 != yCell2)
            zShiftBool = (zCell1 != zCell2)
            return (xShiftBool, yShiftBool, zShiftBool)
        
        
        #Now repeat process for inter-cellular bonds
        for bond in interCellBonds.list:
            pos1 = bond.pos1
            pos2 = bond.pos2
            print "processing intercell bond: ", pos1, " , ", pos2
            pos1Index = -1
            pos2Index = -1
            for i in range(len(allAtoms)):
                currentPos = allAtoms[i].pos
                #if currentPos == pos1:
                if currentPos[0] == pos1[0] and currentPos[1] == pos1[1] and currentPos[2] == pos1[2]:
                    pos1Index = i
                    break
                
            for i in range(len(allAtoms)):
                currentPos = allAtoms[i].pos  
                #if currentPos == pos2:
                if currentPos[0] == pos2[0] and currentPos[1] == pos2[1] and currentPos[2] == pos2[2]:
                    pos2Index = i
                    break
            
            if pos1Index < 0 or pos2Index < 0:
                print "Atom list does not contain all atoms!"
                if pos1Index < 0:
                    print pos1, " missing"
                if pos2Index < 0:
                    print pos2, " missing"
                raise Exception("Export Failed")
            else:
                direction = bondDirection(pos1, pos2)
                allAtoms[pos1Index].addInterCellInteraction(pos2Index, bond.jMatrix, direction)
                allAtoms[pos2Index].addInterCellInteraction(pos1Index, bond.jMatrix, direction)
                print "adding interaction between atom ", pos1Index, " and atom ", pos2Index
        
        
        
#        timer.printTime()
        print"translating..."
        
        def validBond(index1, index2, direction):
            """1/20/11 I beleive the point of this function was to find the interactions that go off
            the edge of the translated lattice and therefore are invalid.  When the larger lattice is
            made, it is made by translating the cutoff cell.  So let's say there is an interaction
            from cell one, to it's neighbo rin the x direction in cell2.  Then we translate in the x
            direction and get cell three, but then we have gone as far as we need to in the x direction
            (made by a for loop) and cell 4 is put next to cell one again, with y coordinate incremented.
            The interaction between cell 3 and cell 4 is invalid becuase they are not next to each other
            in the x direction.  That is the problem, but this method is not working correctly.
            To work correctly it needs to check each direction that the bond crosses a cell boundary in 
            (the direction argument) to see if it would be crossing the edge of the whole translated lattice.
            
            Right now it checks if the bond is in a different translated row of cells, but the problem
            is that legitimate interactions can do this."""
            #print "?valid bond: ", allAtoms[index1].pos, " , ", allAtoms[index2].pos, direction
            cell1 = index1/numAtomsPerCell
            cell2 = index2/numAtomsPerCell
            #Find the coordinates of the cell in units of interaction cells
            posInX1 = int(cell1/(size*size))
            posInX2 = int(cell1/(size*size))
            leftover1 = cell1%(size*size)
            leftover2 = cell2%(size*size)
            posInY1 = int(leftover1/size)
            posInY2 = int(leftover2/size)
            posInZ1 = leftover1%size
            posInZ2 = leftover2%size
            
            #Now, a valid interaction can cross an interaction cell boundary in any direction,
            #but it has a maximum length of one interaction cell.  However, I have made the minimum
            #size of this larger translated lattice equal to 3*3*3 interaction cells.  Therefore,
            #when we hit an edge and get in invalid interaction, the cells will be at least 2
            #interaction cells apart in the direction of the interaction.
            if(direction[0]):
                if numpy.abs(posInX1 - posInX2)>1:
                    #print "false"
                    return False
            if(direction[1]):
                if numpy.abs(posInY1 - posInY2)>1:
                    #print "false"
                    return False
            if(direction[2]):
                if numpy.abs(posInZ1 - posInZ2)>1:
                    #print "false"
                    return False
            print #"true"
            return True

            #Old (incorrect) method:
            if 0:
                print "?valid bond: ", allAtoms[index1].pos, " , ", allAtoms[index2].pos, direction
                cell1 = index1/numAtomsPerCell
                cell2 = index2/numAtomsPerCell
                zRow1 = cell1/size#this relies on the list being created in the nested for loop that was used, z within y within x
                zRow2 = cell2/size
                if(zRow1 != zRow2 and direction[2]):
                    print "false"
                    return False
                xLayer1 = cell1/(size*size)
                xLayer2 = cell2/(size*size)
                if(xLayer1 != xLayer2 and direction[1]):
                    print "false"
                    return False
                #shouldn't have to check z, because if it's not valid in z direction, it would be off the list (>len(allAtoms))
                print "true"
                return True
   
        
        #translate bonds contained within the CutoffCell
        for i in range(len(allAtoms)- numAtomsPerCell):
            newIndex = i+numAtomsPerCell
            for interaction in allAtoms[i].interactions:
                newInteraction = interaction[0] + numAtomsPerCell
                if newInteraction < len(allAtoms):#Should always be the case now
                    allAtoms[newIndex].addInteraction(newInteraction, interaction[1])
                else:#for testing
                    print "\n\ncellbonds contains inter-cutoff cell bond!\n\n"
                    raise Exception("cellbonds contains inter-cutoff cell bond!")
        
            
            
        
        #translate bonds between Cutoff cells
        #size^3 * numAtomsPerCell = len(allAtoms)
        
        #This method iterates through the whole list of atoms.  Each time it encounters
        #an interaction it translates it to all later corresponding indices.  This was
        #a necessary change from the method above, because with the method above, a bond
        #would stop propagating as soon as it encountered one invalid location (an edge).
        #This new method, however, will re-copy interactions that were copied early on 
        #in the main loop, but are encountered later again in the main loop.  This could
        #become slow with large lists.  An alternate method would be to create a copy of
        #the list and copy only from the original to the copy, which eliminates the need
        #for checking repeats and ensures that each interaction is only propagated once.
        cubeSize = size*size*size
        for cell in range(cubeSize):
            for i in range(numAtomsPerCell):
                atomIndex = cell*numAtomsPerCell + i
                for interaction in allAtoms[atomIndex].interCellInteractions:
                    for n in range(1, cubeSize - cell):
                        displacement =  numAtomsPerCell*n
                        if interaction[0] + displacement < len(allAtoms):
                            if validBond(atomIndex + displacement, interaction[0] + displacement, interaction[2]):
                                #Checks for duplicates
                                allAtoms[atomIndex + displacement].addInterCellInteraction(interaction[0] + displacement, interaction[1], interaction[2])
                        
                        
#                    newInteraction = interaction[0] + numAtomsPerCell
#                    if newInteraction < len(allAtoms):
#                        allAtoms[i+numAtomsPerCell].addInterCellInteraction(newInteraction, interaction[1])
    
            
        
#        print "done translating, checking list"
#        timer.printTime()
        
        #Check for reapeats in finalBond List just for testing
#        def isRepeat(finalBondList):
#            for i in range(0, len(finalBondList)):
#                for j in range(i + 1, len(finalBondList)):
#                    if finalBondList[i].sameBond(finalBondList[j]):
#                        return True
#            return False
#        
#        if isRepeat(finalBondList):
#            print "There is a repeat!"
#        else:
#            print "NO repeats!"
#            
#        timer.printTime()
        
        
        #Check the simple atom list
        def atomBalanced(atomIndex):
            atom = allAtoms[atomIndex]
            for otherAtomIndex in range(len(allAtoms)):
                otherAtom = allAtoms[otherAtomIndex]
                if atomInteractsWithAtom(atomIndex, otherAtom):
                    if atomInteractsWithAtom(otherAtomIndex, atom):
                        return True
                    else:
                        return False
            return False
        
        def atomInteractsWithAtom(atomIndex, otherAtom):
            for interaction in otherAtom.interactions:
                if atomIndex == interaction[0]:
                    return True
            return False
        
        
#        for atomIndex in range(len(allAtoms)):
#            if not atomBalanced(atomIndex):
#                print "Not Balanced!!!"
#                break
#        else:
#            print "Balanced!"

        return matrices, allAtoms 
      
    def loadSpinFile(self, fileName):
        """Loads the file output from the simulated annealing that lists the spin
        for each atom."""

        def getSpin(position, lines):
            """Searches lines(a list of the lines found in the spin file) for
            a certain atom and returns the spin associated with it, or None if
            the atom is not in the list of lines."""
            x = str(position[0])
            y = str(position[1])
            z = str(position[2])
            
            spin = None
            
            for line in lines:
                if not line.startswith('#'):
                    values = line.split()
                    if x == values[1] and y == values[2] and z == values[3]:
                        spin = (float(values[4]), float(values[5]), float(values[6]))
                        break
            
            return spin
        
        
        #open the file and creat a list of lines
        file = open(fileName, 'r')
        lines = file.readlines()
        file.close()
        
        for atom in self.getCutoffCell().getAllAtoms():
            newSpin = getSpin(atom.getPosition(), lines)
            if newSpin != None:
                atom.setSpin(newSpin)
                
        #for test purposes
#        print "test:"
#        for atom in self.getCutoffCell().getAllAtoms():
#            print atom.getSpin()
        
        send(signal = "Model Change", sender = "Session")
        
                
            


class Timer():
    """Used to measure the time it takes for various things to compute."""
    def __init__(self):
        self.initialTime = time.clock()
    
    def printTime(self):
        seconds = time.clock() - self.initialTime
        minutes = int(seconds)/60
        seconds -= (minutes*60)
        hours = minutes/60
        minutes -= hours*60
        print "Time:", hours, "hours", minutes, "minutes", seconds, "seconds" 
 
 
    
class atomTable(wx.grid.PyGridTableBase):
    """Contains the information the user enters about atoms. The GUI displays
    what is in this table."""
    def __init__(self):
        wx.grid.PyGridTableBase.__init__(self)
        self.colLabels = ['Symbol', "Atomic Number\n  (optional)",'Valence', '   x   ', '   y   ','   z   ', '  Dx  ', '  Dy  ', '  Dz  ', 'Spin Magnitude']
        self.rowLabels=['Atom 1']
        
        self.data = [
                     ['','','','','', '', 0.0, 0.0, 0.0, '']#Row 1
                     ]
    
    def GetNumberRows(self):
        return len(self.data)
    def AppendRows(self, num):
        for i in range(num):
            self.AppendRow()
        return True
    def GetNumberCols(self):
        return len(self.colLabels)
    def GetColLabelValue(self, col):
        return self.colLabels[col]
    def GetRowLabelValue(self, row):
        return self.rowLabels[row]
    def IsEmptyCell(self, row, col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True
    # Get/Set values in the table.  The Python version of these
    # methods can handle any data-type, (as long as the Editor and
    # Renderer understands the type too,) not just strings as in the
    # C++ version.
    def GetValue(self, row, col):
        try:
            return self.data[row][col]
        except IndexError:
            return ''
    def SetValue(self, row, col, value):
        try:
            self.data[row][col] = value
        except IndexError:
            # add a new row
            self.AppendRow()
            self.data[row][col]=value
        return
    
    def AppendRow(self):
            self.data.append([''] * self.GetNumberCols())
            self.rowLabels.append('Atom ' + str(self.GetNumberRows()))
            self.SetValue(self.GetNumberRows() - 1, 6, 0.0)
            self.SetValue(self.GetNumberRows() - 1, 7, 0.0)
            self.SetValue(self.GetNumberRows() - 1, 8, 0.0)

            # tell the grid we've added a row
            msg = wx.grid.GridTableMessage(self,            # The table
                    wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED, # what we did to it
                    1                                       # how many
                    )
            self.GetView().ProcessTableMessage(msg)
            return True

    def DeleteRows(self,pos=0,numRows=1):
        if numRows>=0 and numRows<=self.GetNumberRows():
            del self.data[pos:pos+numRows]
            msg = wx.grid.GridTableMessage(self,            # The table
            wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED, # what we did to it
            pos,numRows                                     # how many
            )
            #msg = wx.grid.GridTableMessage(self, 0, numRows)
            self.GetView().ProcessTableMessage(msg)
            
            self.UpdateValues()
            return True
        else:
            return False
        
    def UpdateValues( self ):
            """Update all displayed values"""
            msg =wx.grid.GridTableMessage(self, wx.grid.GRIDTABLE_REQUEST_VIEW_GET_VALUES)
            self.GetView().ProcessTableMessage(msg)



class bondTable(wx.grid.PyGridTableBase):
    """Contains the information the user enters about interactions. The GUI displays
    what is in this table."""
    def __init__(self, parameter_manager):
        wx.grid.PyGridTableBase.__init__(self)
        self.paramManager = parameter_manager
        self.colLabels = ['Atom1 Number', '  Na  ','  Nb  ', '  Nc  ', 'Atom2 Number', '  Na  ','  Nb  ', '  Nc  ', '  Jij Matrix  ', 'On']
        self.rowLabels=['Bond 1']
        
        self.data = [
                     ['','','','','','','','','','']#Row 1
                     ]
        
        #Set defualt Jij to 0's  JParam() defaults to 0
        self.SetValue(0, 8, numpy.array([[JParam(self.paramManager),JParam(self.paramManager),JParam(self.paramManager)],
                                         [JParam(self.paramManager),JParam(self.paramManager),JParam(self.paramManager)],
                                         [JParam(self.paramManager),JParam(self.paramManager),JParam(self.paramManager)]]))
        
        #This control what is returned by GetValue for the purposes of displaying.
        #If it is true, the value or range of values of each parameter will be shown.
        #If it is false, the name ('p' + index) will be returned by GetValue
        self.valueView = False
        
        #Each time a Jmatrix is added or removed to the list, the parameter panel needs
        #to be notified.
        self.param_panel = None
    
    def __deepcopy__(self):
        print self.data[0][8]
        for param in self.paramManager.parameters:
            print param.fit
        newParamManager = ParamManager()
        newTable = bondTable(newParamManager)#This will add 9 new JParam objects to the manager
        newParamManager = ParamManager()
        newTable.paramManager = newParamManager
        newTable.data = []
        numParams = len(self.paramManager.parameters)
        #Add all the new Parameters to the new manager
        for i in range(numParams):
            oldParam = self.paramManager.parameters[i]
            JParam(newParamManager, fit = oldParam.fit, value = oldParam.value, min = oldParam.min, max = oldParam.max, default = oldParam.default)
        #Now tie the appropriate parameters
        print "oldParams:\n"
        for param in self.paramManager.parameters:
            print param.group
        print "\n\nnewParams:", len(newParamManager.parameters), "\n\n"
        for index in range(numParams):
            newParam = newParamManager.parameters[index]
            oldParam = self.paramManager.parameters[index]
            for otherParam in oldParam.tied:#These are indices and should match between new and old lists
                newParam.tieTo(otherParam)
        print "\n\nnewParams:", len(newParamManager.parameters), "\n\n"
        #duplicate the data list
        for i in range(len(self.data)):
            newTable.data.append(['','','','','','','','','',''])
            newTable.SetValue(i,1, self.data[i][1])
            newTable.SetValue(i,2, self.data[i][2])
            newTable.SetValue(i,3, self.data[i][3])
            newTable.SetValue(i,4, self.data[i][4])
            newTable.SetValue(i,5, self.data[i][5])
            newTable.SetValue(i,6, self.data[i][6])
            newTable.SetValue(i,7, self.data[i][7])
            newTable.SetValue(i,9, self.data[i][9])
            #Create Matrices with the new parameters in the appropriate positions
            j11 = newParamManager.parameters[self.data[i][8][0][0].GetIndex()]
            j12 = newParamManager.parameters[self.data[i][8][0][1].GetIndex()]
            j13 = newParamManager.parameters[self.data[i][8][0][2].GetIndex()]
            j21 = newParamManager.parameters[self.data[i][8][1][0].GetIndex()]
            print "\n\nnewParams:", len(newParamManager.parameters), "\n\n"
            j22 = newParamManager.parameters[self.data[i][8][1][1].GetIndex()]
            j23 = newParamManager.parameters[self.data[i][8][1][2].GetIndex()]
            j31 = newParamManager.parameters[self.data[i][8][2][0].GetIndex()]
            j32 = newParamManager.parameters[self.data[i][8][2][1].GetIndex()]
            j33 = newParamManager.parameters[self.data[i][8][2][2].GetIndex()]
            print "\n\nnewParams:", len(newParamManager.parameters), "\n\n"
            newTable.SetValue(i,8,numpy.array([[j11,j12,j13],
                                               [j21,j22,j23],
                                               [j31,j32,j33]]))
            
        print "\n\nnewParams:", len(newParamManager.parameters), "\n\n"
        for param in newParamManager.parameters:
            print param.group
        print "original data:\n", self.data
        print "\n\nnewdata:\n", newTable.data
        print "old groups: ", self.paramManager.groups
        print "new groups: ", newParamManager.groups
        return newTable
    
    def GetNumberRows(self):
        return len(self.data)
    def AppendRows(self, num):
        for i in range(num):
            self.AppendRow()
        return True
    def GetNumberCols(self):
        return len(self.colLabels)
    def GetColLabelValue(self, col):
        return self.colLabels[col]
    def GetRowLabelValue(self, row):
        return self.rowLabels[row]
    def IsEmptyCell(self, row, col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True
        except ValueError: #Jij Matrix
            return False
    # Get/Set values in the table.  The Python version of these
    # methods can handle any data-type, (as long as the Editor and
    # Renderer understands the type too,) not just strings as in the
    # C++ version.
    def GetValue(self, row, col):
        """Returns a String representation of the value in the given cell."""
        if self.valueView or col != 8:
            try:
                return str(self.data[row][col])#This must be a string or it cuases problems one the display side
            except IndexError:
                return ''
        else: #display names of paramters
            jStr = ''
            for i in range(3):
                jStr += '['
                for j in range(3):
                    jStr += self.data[row][col][i][j].getName()
                    if j< 2:
                        jStr += "   "
                jStr += ']'
                if i<2:
                    jStr += '\n'
            return jStr
                    
                    
        
    def GetActualValue(self, row, col):
        """Returns the actual value in the cell, not just the string representation
        like GetValue"""
        try:
            return self.data[row][col]
        except IndexError:
            return None
    def SetValue(self, row, col, value):
        try:
            self.data[row][col] = value
        except IndexError:
            # add a new row
            self.AppendRow()
            self.data[row][col]=value
        return
    
 #   def AppendRows(self, numRows):
 #       for i in range(numRows):
 #           self.AppendRow()
 #       return True
    
    def AppendRow(self):
        row = [''] * (self.GetNumberCols())
        self.data.append(row)
        self.rowLabels.append('Bond ' + str(self.GetNumberRows()))
        
        #Set defualt Jij to 0's
        #SetValue does not work for some reason:
        #self.SetValue(self.GetNumberRows()-1, 8, numpy.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]))
        self.data[self.GetNumberRows()-1][8] = numpy.array([[JParam(self.paramManager),JParam(self.paramManager),JParam(self.paramManager)],
                                                            [JParam(self.paramManager),JParam(self.paramManager),JParam(self.paramManager)],
                                                            [JParam(self.paramManager),JParam(self.paramManager),JParam(self.paramManager)]])

        # tell the grid we've added a row
        if self.GetView() != None:
            msg = wx.grid.GridTableMessage(self,            # The table
                    wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED, # what we did to it
                    1                                       # how many
                    )
            self.GetView().ProcessTableMessage(msg)
        
        if self.param_panel:
            self.param_panel.AddRow(self.GetNumberRows()-1)
        
        return True
    
    def DeleteRows(self,pos=0,numRows=1):
        if numRows>=0 and numRows<=self.GetNumberRows():
            #Update teh parameter panel
            self.param_panel.RemoveRows(pos, numRows)
            
            #remove the parameter objects from the manager
            for i in range(pos, pos + numRows):
                for j in range(3):
                    for k in range(3):
                        self.paramManager.removeParam(self.data[i][8][j][k])
            del self.data[pos:pos+numRows]
            msg = wx.grid.GridTableMessage(self,            # The table
            wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED, # what we did to it
            pos,numRows                                     # how many
            )
            #msg = wx.grid.GridTableMessage(self, 0, numRows)
            self.GetView().ProcessTableMessage(msg)
            
            self.UpdateValues()
            return True
        else:
            return False
        
    def UpdateValues( self ):
        """Update all displayed values"""
        msg =wx.grid.GridTableMessage(self, wx.grid.GRIDTABLE_REQUEST_VIEW_GET_VALUES)
        self.GetView().ProcessTableMessage(msg)
        
    def GetJMatrix(self, row):
        return self.GetActualValue(row, 8)
    
    def AddParamPanel(self, panel):
        self.param_panel = panel
        for i in range(self.GetNumberRows()):
            self.param_panel.AddRow(i)
        
