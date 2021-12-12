#  prebiotic synthesis planetesimal

[![DOI](https://zenodo.org/badge/429479674.svg)](https://zenodo.org/badge/latestdoi/429479674)

This code allows calculating the abundances of prebiotic molecules inside planetesimals in the early solar system, which could provide a potentially habitable world with key ingredients for the formation of life. Simulation in **C++** (which is using the proprietary **FORTRAN** library **ChemApp**, not included in repo), **Python3**, **R**, and a **Makefile** for easy building.

[Klaus Paschek](https://www.mpia.de/institute/staff/113334) [![ORCID](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-2603-4236), Max Planck Institute for Astronomy, Heidelberg, Germany

## Table of Contents
1. [ChemApp](#chemapp)
2. [Installation](#installation)
3. [How to run](#how-to-run)
4. [Target molecules](#target-molecules)
5. [Planetesimal temperatures](#planetesimal-temperatures)
6. [Resulting plots](#resulting-plots)


## ChemApp
This Code needs the proprietary **FORTRAN** library [ChemApp](https://gtt-technologies.de/software/chemapp/) to function. It is not provided in this repository. This software is distributed by _GTT-Technologies_ and provided by this company as binaries. A free light version with limited capabilities is available, the full functionality allowing for more chemical substances and elements is available as a paid version. For the chemical reactions provided here, the paid full version is necessary.

## Installation
After cloning this repository, the **ChemApp** binaries are needed. In the **Makefile** the directory containing the binaries has to be specified using the variable `LIBSDIR`. Open the **Makefile** and change the variable accordingly, e.g.:
```
LIBSDIR = ./chemapp-v<your_version>/
```
Further, the [gcc](https://gcc.gnu.org/) compiler is needed for mixed compiling of the **FORTRAN** libraries with the **C++** code. This is done using the flag `-lgfortran` and already set up in the **Makefile**. At least a _C++11_-standard compatible version of **gcc** is needed. **gcc** is usually preinstalled on Linux or can be installed (including **make**) using the command:
```
sudo apt install build-essential
```
For the use of the [CHNOSZ](chnosz.net) database one needs to install **R**:
```
sudo apt-get install r-base r-base-dev
```
To install **CHNOSZ** one needs to start **R** in the terminal and execute the following commands:
```sh
you@your_machine$ R
> install.packages("CHNOSZ")
> q()
Save workspace image? [y/n/c]: n
```
Further, on needs **Python3** (usually already preinstalled on Linux), e.g.:
```
sudo apt install python3.5
```
For the **Python3** code to work one needs the following libraries:

- `pybind11`
- `numpy`
- `scipy`
- `matplotlib`
- `pandas`
- `rpy2`
- `math`
- `datetime`
- `calender`
- `datetime`
- `calendar`
- `pathlib`
- `sys`
- `csv`
- `tqdm`
- `contextlib`
- `warnings`

If not installed yet, they can be installed using **pip3**, e.g.:
```
pip3 install pybind11
```

## How to run
First, one has to compile the **C++** code together with **ChemApp** into a **Cython** module to be usable in **Python3**. This can be done easily using the **Makefile** by simply typing in the directory with the cloned directory:
```
make
```
This will automatically build the software using **gcc** and **pybind11**.
To run a reaction one has to run now the **Python3** script:
```
python3 prebiotic_synthesis.py <target_molecule>[/<target_molecule_2>] [<reaction_no.> OR all] [<pressure_in_bar>]
```
E.g.:
```
python3 prebiotic_synthesis.py adenine 1 100
```
Several _target molecules_ can be given separated by `/` to start the abundance calculation for all of them simultaneously:
```
python3 prebiotic_synthesis.py adenine/ribose 1 100
```
The default _pressure_ is set to _100 bar_ and can be used by just giving no explicit value:
```
python3 prebiotic_synthesis.py adenine 1
```
Some _reaction numbers_ are already set up in the script and one can run all these reactions at once by using either `all`:
```
python3 prebiotic_synthesis.py adenine/ribose/cytosine all 100
```
Or by just neglecting the _reaction number_ completely (be aware that now one can't give an explicit value for the _pressure_ and the default _100 bar_ will be used):
```
python3 prebiotic_synthesis.py adenine/ribose/cytosine
```

### Target molecules
The _target molecule(s)_ has/have to be set up (each) as a **csv** file in the subdirectory `./reaction_info/` with the filename as `<target_molecule>_<reaction no.>.csv`, e.g., `adenine_1.csv`.
In these **csv** files the reactants and products have to be listed together with their phase and initial concentration (normalized to water), e.g.:
```sh
H2O,product
CO,gas,0.00000175
H2,gas,0.00000175
NH3,gas,0.007
H2O,liq,0.0
adenine,aq,0.0
adenine,cr,0.0
```
Molecule names have to follow the naming convention of **CHNOSZ**.
The first line defines the role of water and the second argument can be either `reactant`, `product`, or `solvent` if water does not take part in the reaction.
The phase of the molecule can be `gas`, `liq` (for water only), `aq` (aqueous --- solved in water), or `cr` (crystal --- solid).

### Planetesimal temperatures
The files giving the thermodynamic temperature profiles over time and for several radii inside the planetesimals have to be given in the subdirectory `./temps_input/` as **csv** files.
In these files, the first row has to define the radii in meters inside the planetesimal, and the first column the time after formation in years. The temperatures have to be given in kelvins.
| time | radius 1 | radius 2 | ... |
| ------ | ------ | ------ | ------ |
| time 1 | temperature time 1 radius 1 | temperature time 1 radius 2 | ... |
| time 2 | temperature time 2 radius 1 | temperature time 2 radius 2 | ... |
| ... | ... | ... | ... |

Example:
| 0 | 50 | 1050 | ... |
| ------ | ------ | ------ | ------ |
| 20 | 161.3 | 161.3 | ... |
| 30 | 161.5 | 161.4 | ... |

Example **csv**:
```sh
0,50,1050,...
20,161.3,161.3,...
30,161.5,161.4,...
```

## Resulting plots
These temperature files are read in and the calculations are starting. The in the **Python3** script provided plotting routines will then give the resulting plots in the subdirectory `./results/` as **pdf** files.
