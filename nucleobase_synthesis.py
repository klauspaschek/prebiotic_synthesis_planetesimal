import numpy as np
import matplotlib.pyplot as plt
import pathlib
import sys
import csv

# Custom python libraries
import Water
import GibbsEnergy

# pybind11 C++/FORTRAN library
import ChemApp



####
# Command line arguments
####

# Check if command line arguments are provided
if len(sys.argv) != 4:
    errorStr =  'Command line argument(s) missing. When running '
    errorStr += 'this program, you have to provide 3 command line '
    errorStr += 'arguments:\npython3 fit_gibbs_energy.py <Target '
    errorStr += 'nucleobase> <Reaction no.> <Pressure in bar>'
    raise ValueError(errorStr);

# Assing command line arguments for later use
targetNucleobase =       sys.argv[1]
reactionNo       =       sys.argv[2]
pressure         = float(sys.argv[3])


####
# Obtain information about reaction from input file
####

# Read input file containing information about reactants involved in reaction
# in './reaction_info' subdirectory
# Directory and file path
inputDir      = pathlib.Path('.') / 'reaction_info'
inputFileName = targetNucleobase + '_' + reactionNo + '.csv'
inputPath     = inputDir / inputFileName
if not inputDir.is_dir():
    inputDir.mkdir()
    inputPath.open()
    errorStr =  'Directory \'./' + str(inputDir) + '\' for input file '
    errorStr += 'containing information about reactants involved in reaction '
    errorStr += 'not found'
    raise NotADirectoryError(errorStr)
if not inputPath.is_file():
    inputPath.open()
    errorStr =  'File \'./' + str(inputPath) + '\' '
    errorStr += 'containing information about reactants involved in reaction '
    errorStr += 'not found'
    raise FileNotFoundError(errorStr)

# Read
reactants = []
with inputPath.open() as f:
    readReactants = csv.reader(f, delimiter = ',')
    for row in readReactants:
        reactants.append(row)

# Find role of water in reaction
waterRole = reactants[0][1]


####
# Create thermochemical input file for ChemApp
####

# Calculate temperature of boiling point of water as upper bound for data
# acquisition in R library 'CHNOSZ'
tempMax = Water.boiling_point(pressure)
# Calculate Gibbs energy coefficients for reactants read from input file using
# data from R library 'CHNOSZ'
coeffs = [{}, {}, {}]
for i in range(1, len(reactants)):
    coeffs = GibbsEnergy.fit(reactants[i][0], reactants[i][1],
                             float(reactants[i][2]), coeffs, pressure, tempMax)

# Write input file for ChemApp and obtain initial concentration in correct
# order
initConcs = GibbsEnergy.ChemApp_file(coeffs, targetNucleobase, reactionNo,
                                     pressure, tempMax)


####
# ChemApp
####

# Error variable to check execution of FORTRAN routines in ChemApp
err = int(0)

# Startup ChemApp
ChemApp.start(err)

# Path of thermochemical input file for ChemApp
fileDir = pathlib.Path('.') / 'input_ChemApp'
if not fileDir.is_dir():
    errorStr =  'Directory \'' + str(fileDir) + '\'necessary for ChemApp '
    errorStr += 'to read input file from not found'
    raise NotADirectoryError(errorStr)
filePath = fileDir / (targetNucleobase + '_' + str(reactionNo) + '_' +
                      str(int(pressure)) + 'bar.dat')
if not filePath.is_file():
    errorStr =  'File \'./' + str(filePath) + '\' '
    errorStr += 'containing thermochemical information about reactants '
    errorStr += 'involved in reaction for ChemApp to read from not found'
    raise FileNotFoundError(errorStr)

# Read thermochemical input file into ChemApp
ChemApp.read_data(str(filePath), err)

# Directory for ChemApp to write log files to
logDir = pathlib.Path('.') / 'logs'
logDir.mkdir(exist_ok = True)

# Open log files
ChemApp.open_file(str(logDir), err);

# Do iteration over different temperatures
temps = np.linspace(273.15, tempMax)
amounts = np.empty(len(temps))
for i in range(len(temps)):
    amounts[i] = ChemApp.amount(targetNucleobase, waterRole, initConcs,
                                temps[i], pressure, str(logDir), err)

# Close log files
ChemApp.close_file(err)

# Save amounts in .csv file
amountsFile = pathlib.Path('.') / (targetNucleobase + '_' + str(reactionNo) +
                                   '_' + str(int(pressure)) +
                                   'bar_amounts.csv')
with amountsFile.open(mode = 'w') as f:
    # Write header
    f.write('Temperature[K],Amount[mol ' + targetNucleobase + '/mol H2O]\n')

    for i in range(len(temps)):
        f.write('{:.2f},{:.5e}\n'.format(temps[i], amounts[i]))

# Print amounts to console
for i in range(len(temps)):
    print('{:2.0f}   {:.2f}  {:.5e}'.format(i, temps[i], amounts[i]))

# Plot
fig, ax = plt.subplots()

ax.plot(temps, amounts * 1e9)
ax.set_yscale('log')
ax.set_ylim(1, np.max(amounts * 1e9) * 2)

plt.savefig(targetNucleobase + '_' + str(reactionNo) + '_' +
            str(int(pressure)) + 'bar_amounts.png')
