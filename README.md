# NLO EW and QCD corrections to p p -> e+ ve mu+ mu- A at the LHC



## About



## Prerequisites

* ```gfortran``` : a GNU Fortran compiler.
* ```LHAPDF6``` : a C++ library interfacing with the parton distribution functions (PDFs). ([Installation instructions](https://lhapdf.hepforge.org/))
* ```Recola-Collier``` : a Fortran library providing the QCD/electroweak virtual one-loop amplitude. ([Installation instructions](https://recola.gitlab.io/recola2/installation.html) We used verison 1 of the package when developing the framework, but version 2 should also work and probably should be recommended.)
* ```HELAS``` : a Fortran library needed to calculate the dipole subtraction terms.



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

When the program is used for a publication, please cite:

S. Actis, A. Denner, L. Hofer, A. Scharf, and S. Uccirati, *Recursive generation of one-loop amplitudes in the Standard Model*, JHEP, 04:037, 2013, [arXiv:1211.6316](https://arxiv.org/abs/1211.6316).

S. Actis, A. Denner, L. Hofer, J.N. Lang, A. Scharf, and S. Uccirati, *RECOLA: REcursive Computation of One-Loop Amplitudes*, Comput. Phys. Commun., 214:140–173, 2017, [arXiv:1605.01090](https://arxiv.org/abs/1605.01090).

A. Denner and S. Dittmaier, *Scalar one-loop 4-point integrals*, Nucl. Phys., B844:199–242, 2011, [arXiv:1005.2076](https://arxiv.org/abs/1005.2076).


## Contact

Please send comments and bug reports to huanfeng@buffalo.edu and let us know if you have questions about utilizing the program.
