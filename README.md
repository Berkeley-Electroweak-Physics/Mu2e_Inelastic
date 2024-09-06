# Mu2e_Inelastic
Mathematica code for analyzing inelastic muon-to-electron conversion on the nuclear target $`^{27}`$Al. Implements the nuclear effective theory developed in Refs. [1,2].

# Contents
This repository contains the Mathematica notebook ```Mu2e_Inelastic_v1.nb```, the Mathematica package file ```Mu2e_inelastic_v1.wl```, and two subfolders. 

```Mu2e_Inelastic_v1.nb``` is the primary script that the user will interact with. Upon execution, it loads the necessary functions from the package file ```Mu2e_Inelastic_v1.wl```. Typical users should not need to modify the package file.

The subfolder ```/Densities/``` contains one-body density matrices that govern the relevant transition amplitudes in Al27. There are 12 total density files, representing 4 transitions (proceeding to the ground and first 3 excited states) and 3 different shell-model interactions (Brown-Wildenthal, USDA, and USDB), see Ref. [2] for details. The file names indicate the interaction and transition. For example, ```Al27_usda_0_2.txt``` contains the density matrices for transitions from the ground state (5/2+, 0.0 MeV) to the second excited state, (3/2+, 1.015 MeV).

The subfolder ```/Mu2e_data/``` contains the simulated electron spectra computed by the Mu2e collaboration [3]. The most relevant files are ```mu2e_response_CE.txt``` and ```mu2e_response_DIO.txt```, which respectively contain the anticipated elastic conversion-electron (CE) signal and the primary background from muon decay-in-orbit (DIO). All other backgrounds are negligible.

# Setup
After cloning the Mu2e_Inelastic repository, the first line of ```Mu2e_Inelastic_v1.nb``` should be edited so that ```InelasticRoot``` points to the root folder ```/Mu2e_Inelastic-main/``` with subfolders ```/Densities/``` and ```/Mu2e_Data/```

# Execution
The script is executed by running ```Evaluate Notebook```. The user is prompted to enter the single-nucleon NRET low-energy constants (LECs) $\tilde{c}_i^\tau$, where $\tau$ denotes the isospin. The code expects the LECs to be dimensionless with respect to the weak scale, as defined in Eq. (17) of Ref. [2]. The LECs are input one at a time in the format $\\{i, ~\tilde{c}_i^0, ~\tilde{c}_i^1\\}$. The user does not need to enter the values of all 16 NRET LECs. Once all non-zero values have been entered, the user can enter $\\{0, ~0, ~0\\}$ to proceed with the calculation. Any unspecified LECs are assumed to be zero.

# Output
Given the specified LECs, the code calculates the $\mu\rightarrow e$ conversion rate and branching ratio for each of the four transitions of interest. Each calculation is repeated with each of the 3 different shell-model interactions. The results are printed in the notebook. The results for the 3 different effective interactions are then used to obtain an average result and estimated theoretical uncertainty. 

Next, the code uses the calculated average branching ratios along with the elastic conversion-electron spectrum computed by the Mu2e collaboration to generate the expected conversion-electron spectrum corresponding to the user-specified CLFV scenario. Separate spectra are produced for each transition, and the combined spectrum is calculated by combining these results. These 5 spectra are written to separate files in the folder  ```/Mu2e_Inelastic-main/Output/```. The 3 columns of these files correspond to (1) the reconstructed electron momentum in MeV/c, (2) the expected # of counts in each bin of width 50 keV/c, (3) the uncertainity in the expected # of counts in each bin of width 50 keV/c.

# Citation
If you use this software in your work, please cite:

```
@article{Haxton:2024ecp,
    author = "Haxton, W. C. and Rule, Evan",
    title = "{Distinguishing charged lepton flavor violation scenarios with inelastic $\mu\rightarrow e$ conversion}",
    eprint = "2404.17166",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "N3AS-24-014",
    month = "4",
    year = "2024"
}

@article{Mu2e:2022ggl,
    author = "Abdi, F. and others",
    collaboration = "Mu2e",
    title = "{Mu2e Run I Sensitivity Projections for the Neutrinoless $\mu^- \to e^-$ Conversion Search in Aluminum}",
    eprint = "2210.11380",
    archivePrefix = "arXiv",
    primaryClass = "hep-ex",
    reportNumber = "FERMILAB-PUB-22-749-PPD",
    doi = "10.3390/universe9010054",
    journal = "Universe",
    volume = "9",
    number = "1",
    pages = "54",
    year = "2023"
}
```
# References
[1] [Haxton, W. C., and Rule, E., hep-ph/2404.17166](https://arxiv.org/abs/2404.17166)

[2] Haxton, W. C., and Rule, E., hep-ph/2409.xxxxx

[3] [Mu2e Collaboration, Universe 9 (2023) 1, 54](https://www.mdpi.com/2218-1997/9/1/54)
