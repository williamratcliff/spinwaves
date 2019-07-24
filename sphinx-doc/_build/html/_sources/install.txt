Installing Spinwaves
====================

Windows Users
-------------

The easiest way for windows users to install the application is with the latest
windows installer available at:

	http://code.google.com/p/spinwaves/downloads/list

Installing From Source
----------------------

You must have python2.4 or python2.6 installed on your computer.

The latest source code can be downloaded from the svn repository::

	$ svn co http://spinwaves.googlecode.com/svn/trunk

Then from the main source directory run::

	$ python setup.py install

This should resolve all dependencies via easy_install and the program can now be run from the shell with::

	$ spinwaves

Dependencies
------------

	The Spinwaves application requires the following third party pieces of software to run:	

	* Python 2.5-2.6
	* wxPython
	* numpy
	* matplotlib
	* sympy
	* scipy
	* vtk
	* Python-Multiprocessing (Python2.5)

	If you are using Python2.5, a back port of the multiprocessing package
	must be installed.  It can be found here:
		http://code.google.com/p/python-multiprocessing/

	Most Linux distributions will come with a version of Python.  To check
	your version, open python in the shell.  The version will be printed.

	Vtk can be acquired from:

		http://www.vtk.org

	All other dependencies should be handled automatically by setuptools
	when the setup.py script is run.  However, since this has proven to be
	error-prone alternate methods of resolving these dependencies have been
	given below.


Ubuntu/Debian/Any Distribution using Aptitude Users
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	A shell script has been included to automatically download all
	dependencies (other than Python itself) using aptitude.  The script
	will download and install the current versions of numpy, matplotlib,
	sympy, scipy, and vtk.  wxPython2.8 will be installed, although 2.6 will
	work as well.  From the top level source folder, type::

	$ bash dependencies.sh

	If you are using Python2.5, a back port of the multiprocessing package
	must be installed.  It can be found here:
		http://code.google.com/p/python-multiprocessing/

