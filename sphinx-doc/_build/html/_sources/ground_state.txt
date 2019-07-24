Finding the Ground State
========================

To calculate the ground state, you must first export the atom interactions to a file, and then run the simulated annealing algorithm from the data in that file.  To export the interaction information to a file, select “Export for Monte Carlo” from the “Compute” menu.  You will be asked how many times you would like to translate the cell for the simulation.  This is asking the size you would like to make the lattice to run the simulation on.  A larger lattice will reduce edge effects, but take longer to compute.

A simulated Annealing algorithm is used to compute the ground state.  This is a type of monte carlo algorithm, where random spins are chosen and those with the lowest resulting total energy are kept.   Each time a new random spin configuration is created, it is either kept or thrown out with a probability depending on the energy and the current “temperature” that the simulation is at.  The simulation begins at a high temperature where there is a high probability of accepting new spin configurations even if they have a higher resulting energy.  The temperature is then slowly decreased, until the system settles into it's (hopefully) global minimum.  This decreasing temperature behaviour is used to reduce the chances of getting stuck in a local minimum.
 
After you have created the interaction file, select “Run Simulation” from the “Compute” menu.  You will be asked to enter the maximum and minimum temperatures as well as the temperature factor and steps per temperature level.  These parameters can be changed for specific cases.  Here are a few things to keep in mind:

	* Maximum temperature:  This is the temperature that the simulation starts at.  If this is too high, it will only make the simulation take a long time.  However, if it is too low, it is likely to get stuck in a local minimum and not find the true ground state.
	* Minimum Temperature:  If this is too high, the Hamiltonian may not have settled into it's minimum and you may get a random configuration.  The local optimizer which is run after the simulation will then return a local minimum value.  If the minimum temperature is too low, the simulation will just take a long time.
	* Steps per temperature level:  This is how many random spin configurations will be tried at each temperature.  A higher number will be more likely to explore all possible values, but will take longer to run.
	* Temperature Factor:  This is the number that the temperature is multiplied by at each step.  A number closer to 1 will slow the cooling process, and therefore make it more reliable, but will also take longer.

If you wish to see the spins in the 3D model, select “Load Spins from file” from the “Compute” menu and select the file you just created in the simulation.
