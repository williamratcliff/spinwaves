#
#This will compile the shared library needed to run the Monte Carlo simulation
#to find the ground state spins and move it to the appropriate directory.
#Scons must be installed on the computer.
#

scons
mv lib_monteCarlo.so ../spinwaves/MonteCarlo/_monteCarlo.so
