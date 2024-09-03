# Mu2e_Inelastic
Mathematica code for analyzing inelastic muon-to-electron conversion on Al27

# Contents
This repository contains the Mathematica notebook Mu2e_Inelastic.nb and two subfolders. 

The subfolder Mu2e_Inelastic/Densities contains one-body density matrices that govern the relevant transition amplitudes in Al27. There are 12 total density files, representing 4 transitions (proceeding to the ground and first 3 excited states) and 3 different shell model interactions (Brown-Wildenthal, USDA, and USDB). The file names indicate the interaction and transition. For example, Al27_usda_0_2.txt contains the density matrices for transitions from the (5/2+, gs) to the (3/2+, 1.015 MeV) excited state.

The subfolder Mu2e_Inelastic/Mu2e_data contains the simulated electron spectra computed by the Mu2e collaboration. 

# Setup
After cloning the Mu2e_Inelastic repository, the first line of Mu2e_Inelastic.nb should be edited so that InelasticRoot points to the root folder Mu2e_Inelastic with subfolders Mu2e_Inelastic/Densities and Mu2e_Inelastic/Mu2e_Data

# Execution
The script is executed by running "Evaluate Notebook". The user is prompted to enter the single-nucleon NRET low-energy constants. The code expects the LECs to be dimensionless with respect to the weak scale, as defined in Eq. (17) of []. Any omitted LECs are assumed to be zero.

# Citation
If you use this software in your work, please cite


