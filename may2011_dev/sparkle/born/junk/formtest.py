from sample_prep import *
import approximations,view, scatter
from numpy import log,abs, min, max
from pylab import show
import time
testcube = Parallelapiped(4.5e-6,dim = [5.0e4,5.0e4,2000.0])
scene = Scene([testcube])

unit = Unit_Cell(Dxyz = [1.0e5,1.0e5,2500.0],n = [50,50,50], scene = scene)
unit.render()

q_space = Q_space([-.0001,-0.0001,0.00002],[.0001,.0001,0.02],[50,50,50])
lattice = Rectilinear([20,20,1],unit)
beam = Beam(5.0,None,None,0.05,None)


calc = scatter.Calculator(lattice,beam,q_space,unit)

t0 = time.time()
calc.cudaSMBA()
print "time",time.time()-t0

a =calc.results
t0 = time.time()
calc.BA()
print "time",time.time()-t0
b =calc.results
#calc.BA()

view.data_compare(a,b,q_space.minimums,q_space.maximums)
show()