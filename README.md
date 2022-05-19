# Automated Monte Carlo Simulator at NLO Accuracy


This is a Monte Carlo (MC) software capable of calculating the next-to-leading order (NLO) QCD and electroweak (EW) corrections to triboson production with fully leptonic decays, ie. p p -> e+ ve mu+ mu- A, at the CERN Large Hadron Collider (LHC).

In this documentation, we would focus on the ***structure***, the ***installation*** and the ***usage*** of the program. For a high-level walk-through of the physics behind the program, please stop by this [article](https://huanfengc.github.io/projects/prediction-of-triboson-production-at-the-lhc/). For the theoretical and technical details of the program, we refer the users/readers to the publication below and the references therein:

> **H. Cheng** and D. Wackeroth, *NLO electroweak and QCD corrections to the production of a photon with three charged lepton plus missing energy at the LHC*, [Phys. Rev. D **105** 096009 (2022)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.105.096009), arXiv: [2112.12052](https://arxiv.org/abs/2112.12052).


This documentation consists of following sections:

* [Prerequisites](#Prerequisites)
* [One more step before running](#One-more-step-before-running)
* [Structure of directories](#Structure-of-directories)
* [Usage](#Usage)
* [Cite](#Cite)
* [Contact](#Contact)
<a name="headers"/>

We suggest that the first-time users read the first four sections in order.



## Prerequisites

To run the program successfully, the following compiler and libraries have to be pre-intalled:

* ```gfortran``` : a GNU Fortran compiler.
* ```LHAPDF6``` : a C++ library interfacing with the parton distribution functions (PDFs). ([Installation instructions](https://lhapdf.hepforge.org/))
* ```Recola-Collier``` : a Fortran library providing the QCD/electroweak virtual one-loop amplitude. We used verison 1 of the package when developing the program, but version 2 should also work and probably should be recommended. ([Installation instructions](https://recola.gitlab.io/recola2/installation.html))
* ```HELAS``` : a Fortran library needed to calculate the dipole subtraction terms. To install  ```HELAS```, please follow the instructions as follows:
  * TBD
  * TBD
  * TBD


## Structure of directories



## One more step before running

```tables.rcl.f90```



## Usage

Download the repository, go into the directory and compile,

```
cd wza-nlo-ew
make
```

run the executable,
```
./epvemumua.x
```
then type in a random seed (eg. '```0```') and press ```enter```. 

When the program is running, the progress of Monte Carlo integration will be printed on the screen. Once the computation is finished, the results will be saved in ```plots/```. In ```plots/```, ```xs.dat``` contains the total cross section and the uncertainty of Monte Carlo integration. All the other ```*.dat``` files are differential cross sections in which the 1st column is the bin, the 2nd column is the differential cross section, the 3rd column is the uncertainty of Monte Carlo integration. 

The unit of cross section calculated by the program is ```femtobarn (fb)```.



## Cite

**When the program is used for a publication, please cite:**

> H. Cheng and D. Wackeroth, *NLO electroweak and QCD corrections to the production of a photon with three charged lepton plus missing energy at the LHC*, [arXiv:2112.12052](https://arxiv.org/abs/2112.12052).

> S. Actis, A. Denner, L. Hofer, A. Scharf and S. Uccirati, *Recursive generation of one-loop amplitudes in the Standard Model*, JHEP, 04, 037 (2013), [arXiv:1211.6316](https://arxiv.org/abs/1211.6316).

> S. Actis, A. Denner, L. Hofer, J.N. Lang, A. Scharf and S. Uccirati, *RECOLA: REcursive Computation of One-Loop Amplitudes*, Comput. Phys. Commun., 214:140–173 (2017), [arXiv:1605.01090](https://arxiv.org/abs/1605.01090).

> A. Denner and S. Dittmaier, *Scalar one-loop 4-point integrals*, Nucl. Phys., B844:199–242 (2011), [arXiv:1005.2076](https://arxiv.org/abs/1005.2076).

> A. Denner and S. Dittmaier, *Reduction of one loop tensor five point integrals*, Nucl. Phys., B658:175–202 (2003), [arXiv:hep-ph/0212259](https://arxiv.org/abs/hep-ph/0212259).

> A. Denner and S. Dittmaier, *Reduction schemes for one-loop tensor integrals*, Nucl. Phys., B734:62–115 (2006), [arXiv:hep-ph/0509141](https://arxiv.org/abs/hep-ph/0509141).

> A. Denner, S. Dittmaier and L. Hofer, *Collier: a fortran-based Complex One-Loop LIbrary in Extended Regularizations*, Comput. Phys. Commun., 212:220–238 (2017), [arXiv:1604.06792](https://arxiv.org/abs/1604.06792).

> T. Gehrmann and N. Greiner, *Photon Radiation with MadDipole*, JHEP 12, 050 (2010), [arXiv:1011.0321](https://arxiv.org/abs/1011.0321).

> R. Frederix, T. Gehrmann and N. Greiner, *Integrated dipoles with MadDipole in the MadGraph framework*, JHEP 06, 086 (2010), [arXiv:1004.2905](https://arxiv.org/abs/1004.2905).

> R. Frederix, T. Gehrmann and N. Greiner, *Automation of the Dipole Subtraction Method in MadGraph/MadEvent*, JHEP 09, 122 (2008), [arXiv:0808.2128](https://arxiv.org/abs/0808.2128).

> A. Buckley, J. Ferrando, S. Lloyd, K. Nordstrom, B. Page, M. Rufenacht, M. Schonherr and G. Watt, *LHAPDF6: parton density access in the LHC precision era*, Eur. Phys. J. C 75, 132 (2015), [arXiv:1412.7420](https://arxiv.org/abs/1412.7420).



## Contact

Please send comments and bug reports to huanfeng@buffalo.edu and let us know if you have questions about utilizing the program.
