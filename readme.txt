Sep 2022 

This code was written to perform the calculations published in:
"Accurate quantification of lattice temperature dynamics from ultrafast electron diffraction of single-crystal films using dynamical scattering simulations" 
by Daniel B. Durham, Colin Ophus, Khalid M. Siddiqui, Andrew M. Minor, and Daniele Filippetto
Preprint posted to ArXiv in Sep 2022.
Published in Structural Dynamics on Dec 05 2022. https://doi.org/10.1063/4.0000170

This code package includes kinematical, bloch wave, and multislice codes to calculate electron diffraction patterns from flat or rippled single crystal foils for UED, with gold included as an example. Some key parameters accounted for here include: crystal structure, lattice temperature, sample topography. Libraries of diffraction patterns can be calculated and fit to experimental data to retrieve lattice temperture changes.

All are welcome to use/adapt this code and the methods described here and in the article. We simply ask if you do that you please cite the article.

We acknowledge support from STROBE: A National Science Foundation Science and Technology Center and the US Department of Energy, Office of Science, Basic Energy Sciences.

Folder organization:
Kinematical: functions specific to kinematical calculations
BlochWaves: functions specific to bloch wave calculations
Multislice: functions specific to multislice calculations
Computations: scripts for performing example computations, including those used in the 2022 ArXiv preprint
Utilities: functions used commonly across the three calculation types
Visualization: functions for generating commonly called for plots
ExperimentFitting: Scripts for comparing simulations to UED experiment, and UED data examples for gold
MaterialFiles: Functions for generating material simulation cells, currently just has the example for gold. 
