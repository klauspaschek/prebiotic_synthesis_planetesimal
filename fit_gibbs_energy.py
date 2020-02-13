import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import pathlib
import pandas as pd
import rpy2.robjects as ro
import rpy2.robjects.packages
import rpy2.robjects.pandas2ri
import math
import datetime
import calendar
import sys
import csv


####
# Calculate approximate boiling point temperature T of water as function of
# pressure P by using an approximation for the pressure temperature dependence
# P(T) = f(T). By rearranging the formula as P(T) - f(T) = g(T) = 0 and finding
# the root using the Newton-Raphson method this formula is solved for T.
####
# Parameters: pressure    : float
#                 Pressure in unit [bar].
####
# Returns:    temperature : float
#                 Temperature of boiling point of water in Kelvin.
def water_boiling_point(pressure):
    # Parameters of approximating function for temperature pressure dependence
    a = -6094.4642
    b = 21.1249952
    c = -2.7245552e-2
    d = 1.6853396e-5
    e = 2.4575506

    # Rearranged function P(T) - f(T)
    def f(T):
        return pressure - np.exp(a / T + b + c * T + d * T**2 + \
                                 e * np.log(T)) * 1e-5

    # Derivative of f(T)
    def fprime(T):
        return -np.exp(a / T + b + c * T + d * T**2 + e * np.log(T)) * \
               (-a / T**2 + c + 2 * d * T + e / T) * 1e-5

    temperature = scipy.optimize.newton(f, 1000, fprime = fprime)

    return temperature


####
# Obtain data from R database 'CHNOSZ'
####
# Parameters: molecule : str
#                 String containing the name of the wanted molecule as stored
#                 in the R 'CHNOSZ' database.
#             phase    : str
#                 String containing the name of the wanted phase of the
#                 molecule as stored in the R 'CHNOSZ' database.
#             tempMax  : float
#                 Temperature up to which data should be acquired in unit
#                 Kelvin.
####
# Returns:    out      : numpy.2darray
#                 Contains temperatures as first and corresponding Gibbs
#                 energies as second element.
def get_R_data(molecule, phase, tempMax):
    tempMaxStr = '{:.2f}'.format(tempMax)

    # Use rpy2 library to start R session and obtain data
    # Load database in R
    ro.packages.importr('CHNOSZ')

    # Change units to Kelvin and bar
    ro.r('T.units(\'K\')')
    ro.r('P.units(\'bar\')')

    # Compose command to obtain data from R as string
    Rcommand =  'subcrt(info(\'' + molecule + '\', state = \'' + phase
    Rcommand += '\')[[1]], T = seq(273.15, ' + tempMaxStr + ', '
    Rcommand += 'length.out = 200), P = ' + str(pressure) + ')'

    # Get data as R dataframe
    hyd = ro.r(Rcommand)

    # Convert to numpy array (with pandas dataframe as itermediate step)
    ro.pandas2ri.activate()
    pd_hyd = ro.conversion.rpy2py(hyd)

    out = np.transpose(pd_hyd[1][0][['T', 'G']].to_numpy())

    # Change energy unit from [cal] to [kJ]
    out[1] = out[1] * 4.1858 / 1e3

    return out


####
# Fit Gibbs energy coefficients and store them in dictionary
####
# Parameters: molecule : str
#                 String containing the name of the wanted molecule as stored
#                 in the R 'CHNOSZ' database.
#             phase    : str
#                 String containing the name of the wanted phase of the
#                 molecule as stored in the R 'CHNOSZ' database.
#             coeffs   : list of dict
#                 List of dictionaries to add new fitted coefficients. First
#                 element contains dictionary for gas phase, second for aqeous
#                 phase and third for condensed phase. Keys of dictionaries are
#                 the names of the molecules as string and values are
#                 numpy.arrays containing the coefficents. Will be returned.
#             pressure : float
#                 Pressure in unit [bar].
#             plot     : boolean, optional
#                 Whether to fit the data and the fit and store them in
#                 './fit_plots' subdirectory. Default is True.
####
# Returns:    coeffs   : list of dict
#                 List of dictionaries to add new fitted coefficients. First
#                 element contains dictionary for gas phase, second for aqeous
#                 phase and third for condensed phase. Keys of dictionaries are
#                 the names of the molecules as string and values are
#                 numpy.arrays containing the coefficents.
def fit_gibbs_energy(molecule, phase, coeffs, pressure, plot = True):
    # Fit function for Gibbs energy
    def gibbs_energy(T, a, b, c, d, e, f):
        return a + b * T + c * T * np.log(T) + d * T**2 + e * T**3 + f / T

    # Obtain data from R database 'CHNOSZ'
    data = get_R_data(molecule, phase, pressure)

    # Fit gibbs energy coefficients
    popt, pcov = scipy.optimize.curve_fit(gibbs_energy, data[0], data[1])

    # Collect them in list of dict
    if phase == 'gas':
        coeffs[0][molecule] = popt
    elif phase == 'aq':
        coeffs[1][molecule] = popt
    elif phase == 'cr':
        coeffs[2][molecule] = popt
    else:
        errStr  = '\'' + phase + '\' as phase argument not allowed, '
        errStr += 'has to be either \'gas\', \'aq\' or \'cr\'!'
        raise ValueError(errStr)

    # Plot data and fit if input parameter 'plot = True'
    if (plot):
        fig, ax = plt.subplots()

        # Data
        ax.scatter(data[0], data[1], s = 7,
                   label = 'Data from CHNOSZ database')

        # Fit
        fitLabel = 'Fit: %1.3f+%1.3fT+%1.3fTlogT+\n%1.3fT^2+%1.3fT^3+%1.3f/T' \
                   % tuple(popt)
        ax.plot(data[0], gibbs_energy(data[0], *popt), color = 'C1',
                label = fitLabel)

        ax.set_xlabel('T [K]')
        ax.set_ylabel('Gibbs energy of formation [kJ/mol]')
        plt.title(molecule + '(' + phase + ') at pressure of 100 bar')
        plt.legend()

        # Compose file path to store plot to
        fitPlotsDir  = pathlib.Path('.') / 'fit_plots'
        # Create subdirectory if necessary
        fitPlotsDir.mkdir(exist_ok = True)
        plotFileName = molecule + '_' + phase + '_100atm_fit.png'
        plotPath     = fitPlotsDir / plotFileName

        # Save plot
        plt.savefig(str(plotPath))

    return coeffs


# Elementary molecule composition
elementsCompData = [{'H': 1, 'C': 1, 'N': 1},    # HCN
                    {'H': 5, 'C': 5, 'N': 5},    # adenine
                    {'H': 2, 'O': 1},            # water
                    {'H': 3, 'N': 1}]            # NH3
elementsCompIndices = ['HCN', 'adenine', 'H2O', 'NH3']
elementsComp = pd.DataFrame(elementsCompData, index = elementsCompIndices,
                            dtype = float)

# Element masses (unit [u])
elementsMass = {'H':  1.008,
                'C': 12.011,
                'N': 14.007,
                'O': 15.999}

####
# Write ChemApp input file and store in './input_ChemApp/' subdirectory
####
# Parameters: coeffs      : list of dict
#                 List of dictionaries containing fitted Gibbs energy
#                 coefficients. First element has to contain dictionary for gas
#                 phase, second for aqeous phase and third for condensed phase.
#                 Keys of dictionaries have to be the names of the molecules as
#                 string and values are numpy.arrays containing the
#                 coefficents.
#             nucleobase  : str
#                 Name of nucleobase reaction results in.
#             reactionNum : int
#                 Number of reaction.
#             tempMax     : float
#                 Temperature up to which data is valid in unit Kelvin.
####
# Returns:    <no return>
def input_file_ChemApp(coeffs, nucleobase, reactionNum, tempMax):
    # Identify for molecules necessary elements
    elementsSet = set()
    for phase in coeffs:
        for molecule in phase:
            for i in range(len(elementsComp.loc[molecule])):
                if not math.isnan(elementsComp.ix[molecule, i]):
                    elementsSet.add(elementsComp.columns.values[i])
    elementsList = list(elementsSet)

    # Path of file
    fileDir = pathlib.Path('.') / 'input_ChemApp'
    # Create subdirectory if necessary
    fileDir.mkdir(exist_ok = True)
    filePath = fileDir / (nucleobase + '_' + str(reactionNum) + '_100bar' +
                          '.dat')

    with filePath.open(mode = 'w') as f:
        ####
        # Header lines
        ####
        # element1-element2-...
        for element in elementsList[:-1]:
            f.write(element + '-')
        f.write(elementsList[-1])
        today = datetime.datetime.now()
        f.write('    See commentary below!    ')
        f.write('MPIA ' + today.strftime('%b') + ' ' + today.strftime('%Y') +
                '\n')

        # Phase information
        # Number of constituents
        f.write(str(len(elementsList)) + '   ')
        # Numbers of mixing phases
        numMixPhases = 0
        for phase in coeffs:
            numMixPhases += len(phase)
        f.write(str(numMixPhases) + '  ')
        for phase in coeffs[:-1]:
            f.write(str(len(phase)) + '  ')
        # Number of condensed phases
        f.write(' ' + str(len(coeffs[-1])) + '\n')

        # Elements (Constituents)
        # Convert string to FORTRAN compatible format
        def to_FORTRAN_str(string):
            while len(string) < 25:
                string += ' '
            return string

        def write_elements(elements):
            for i in range(len(elements)):
                if i % 3 != 0 or i == 0:
                    f.write(to_FORTRAN_str(elements[i]))
                else:
                    f.write('\n' + to_FORTRAN_str(elements[i]))

        # Write element names
        write_elements(elementsList)
        f.write('\n')
        # Find atomic masses for elements
        elementsMassList = []
        for element in elementsList:
            elementsMassList.append(str(elementsMass[element]))
        # Write element masses (in unit [u])
        write_elements(elementsMassList)
        f.write('\n')

        # Gibbs energy equation info
        # Temperature dependence
        # G(T) = a + b * T + c * T * ln(T) + d * T^2 + e * T^3 + f / T
        # This line defines, that all six coefficients from 'a' - 'f' are used
        f.write('6   1   2   3   4   5   6\n')
        # Pressure dependence
        # This line defines that there is no pressure dependence, so
        # G(P) = a
        f.write('1   1\n')

        ####
        # Mixing Phases
        ####
        for i in range(len(coeffs)):
            if len(coeffs[i]) > 0:
                if i == 0:
                    f.write(to_FORTRAN_str('GAS') + '\n' +
                            to_FORTRAN_str('IDMX') + '\n')
                elif i == 1:
                    f.write(to_FORTRAN_str('AQUEOUS') + '\n' +
                            to_FORTRAN_str('IDMX') + '\n')

                for molecule in coeffs[i]:
                    f.write(to_FORTRAN_str(molecule.upper()) + '\n')
                    f.write('1  1   ')
                    for element in elementsList:
                        f.write(str(np.nan_to_num(elementsComp.loc[molecule,
                                                                   element])) +
                                '   ')
                    f.write('\n')

                    # Expand string with spaces to length of 15
                    def to_15_str(string):
                        while len(string) < 15:
                            string += ' '
                        return string

                    def write_gibbs_coeffs(coeffsList, tempMax):
                        coeffsList = np.append(tempMax, coeffsList)
                        for i in range(len(coeffsList)):
                            if i % 5 != 0 or i == 0:
                                valueStr = '{:.5f}'.format(coeffsList[i])
                                f.write(to_15_str(valueStr))
                            else:
                                valueStr = '{:.5f}'.format(coeffsList[i])
                                f.write('\n' + to_15_str(valueStr))

                    write_gibbs_coeffs(coeffs[i][molecule], tempMax)
                    f.write('\n')


# Check if command line arguments are provided
if len(sys.argv) != 4:
    errorStr =  'Command line argument(s) missing. When running '
    errorStr += 'this program, you have to provide 3 command line '
    errorStr += 'arguments:\npython3 fit_gibbs_energy.py <Target '
    errorStr += 'nucleobase> <Reaction no.> <Pressure in bar>'
    raise ValueError(errorStr);

# Command line arguments
targetNucleobase =       sys.argv[1]
reactionNo       =       sys.argv[2]
pressure         = float(sys.argv[3])

# Calculate temperature of boiling point of water as upper bound for data
# acquisition
tempMax = water_boiling_point(pressure)

# Read input file containing information about reactants involved in reaction
# in './reaction_info' subdirectory
inputDir = pathlib.Path('.') / 'reaction_info'
if not inputDir.is_dir():
    errorStr =  'Directory \'./' + str(inputDir) + '\' for input file '
    errorStr += 'containing information about reactants involved in reaction '
    errorStr += 'not found'
    raise NotADirectoryError(errorStr)
inputFileName = targetNucleobase + '_' + reactionNo + '_reactants.dat'
inputPath     = inputDir / inputFileName
if not inputPath.is_file():
    errorStr =  'File \'./' + str(inputPath) + '\' '
    errorStr += 'containing information about reactants involved in reaction '
    errorStr += 'not found'
    raise FileNotFoundError(errorStr)

reactants = []
with inputPath.open() as f:
    readReactants = csv.reader(f, delimiter = ',')
    for row in readReactants:
        reactants.append([row[0], row[1]])

# Calculate Gibbs energy coefficients for reactants read from input file
coeffs = [{}, {}, {}]
for i in range(1, len(reactants)):
    coeffs = fit_gibbs_energy(reactants[i][0], reactants[i][1],
                              coeffs, tempMax)

#coeffs = fit_gibbs_energy('HCN', 'aq', coeffs, tempMax)
#coeffs = fit_gibbs_energy('adenine', 'aq', coeffs, tempMax)
#coeffs = fit_gibbs_energy('H2O', 'aq', coeffs)
#coeffs = fit_gibbs_energy('adenine', 'cr', coeffs)
#coeffs = fit_gibbs_energy('NH3', 'gas', coeffs)

input_file_ChemApp(coeffs, targetNucleobase, reactionNo, tempMax)
