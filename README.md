# NLO EW and QCD corrections to p p -> e+ ve mu+ mu- A at the LHC


## Prerequisites

* ```gfortran```
* ```LHAPDF``` library
* ```RECOLA``` library
* ```HELAS``` library

## Installation

```
cd wza-nlo-ew
make
```

## Usage

Run the executable,
```
./epvemumua.x
```
then type in a random seed (eg. '```0```') and press ```enter```. 

When the program is running, the progress of Monte Carlo integration will be printed on the screen. Once the computation is finished, the results will be saved in ```plots/``` where ```xs.dat``` contains the total cross section and the uncertainty of Monte Carlo integration, all the other ```*.dat``` files are differential cross sections.
