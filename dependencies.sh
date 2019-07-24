#!/bin/bash

#This shell script will install all dependencies on an Ubuntu or Debian machine.

echo 'installing Numpy'
sudo aptitude install python-numpy
echo 'installing vtk'
sudo aptitude install python-vtk
echo 'installing matplotlib'
sudo aptitude install python-matplotlib
echo 'installing sympy'
sudo aptitude install python-sympy
echo 'installing scipy'
sudo aptitude install python-scipy
echo 'installing wxPython 2.8'
sudo aptitude install python-wxgtk2.8

echo 'If wxPython 2.6 is installed, it will be removed because matplotlib required wxPython2.8 and has errors when two versions are installed.'
sudo aptitude remove python-wxgtk2.6

