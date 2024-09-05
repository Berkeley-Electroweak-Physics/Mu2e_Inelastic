# Mu2e_Inelastic
Mathematica code for analyzing inelastic muon-to-electron conversion on Al27

# Contents
This repository contains the Mathematica notebook ```Mu2e_Inelastic_v1.nb``` and two subfolders. 

The subfolder ```Mu2e_Inelastic/Densities``` contains one-body density matrices that govern the relevant transition amplitudes in Al27. There are 12 total density files, representing 4 transitions (proceeding to the ground and first 3 excited states) and 3 different shell model interactions (Brown-Wildenthal, USDA, and USDB). The file names indicate the interaction and transition. For example, ```Al27_usda_0_2.txt``` contains the density matrices for transitions from the ground state (5/2+, 0.0) to the second excited state, (3/2+, 1.015 MeV).

The subfolder ```Mu2e_Inelastic/Mu2e_data``` contains the simulated electron spectra computed by the Mu2e collaboration. 

# Setup
After cloning the Mu2e_Inelastic repository, the first line of ```Mu2e_Inelastic_v1.nb``` should be edited so that ```InelasticRoot``` points to the root folder ```Mu2e_Inelastic``` with subfolders ```Mu2e_Inelastic/Densities``` and ```Mu2e_Inelastic/Mu2e_Data```

# Execution
The script is executed by running "Evaluate Notebook". The user is prompted to enter the single-nucleon NRET low-energy constants (LECs) $\tilde{c}_i^\tau$, where $\tau$ denotes the isospin. The code expects the LECs to be dimensionless with respect to the weak scale, as defined in Eq. (17) of []. They are input one at a time in the format $\\{i, ~\tilde{c}_i^0, ~\tilde{c}_i^1\\}$. One does not need to enter the values of all 16 NRET LECs. Once all non-zero values have been entered, the user can enter $\\{0, ~0, ~0\\}$ to proceed with the calculation. Any unspecified LECs are assumed to be zero.

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

