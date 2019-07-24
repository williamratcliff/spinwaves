import numpy
import matplotlib.pyplot as plt

a = numpy.loadtxt("times.txt")
plt.plot(a)
plt.show()