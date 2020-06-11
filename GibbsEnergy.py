"""Obtain data for Gibbs energies of molecules from the R library 'CHNOSZ' for
different temperatures T and pressures P, fit the Gibbs coefficients to the
data following
G(T) = a + b*T + c*T*log(T) + d*T^2 + e*T^3 + f/T
and write them to ChemApp compatible data file.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import pandas as pd
import pathlib
import rpy2.robjects as ro
import rpy2.robjects.packages
import rpy2.robjects.pandas2ri
import math
import datetime
import calendar



# Elementary molecule composition
elementsCompData = [{'H': 1, 'C': 1, 'N': 1        }, # HCN
                    {'H': 5, 'C': 5, 'N': 5        }, # adenine
                    {'H': 2,                 'O': 1}, # H2O
                    {'H': 3,         'N': 1        }, # NH3
                    {'H': 2                        }, # H2
                    {        'C': 1,         'O': 1}, # CO
                    {'H': 2, 'C': 1, 'N': 1, 'O': 1}, # CH2NO
                    {'H': 2, 'C': 1,         'O': 1}, # formaldehyde
                    {'H': 4, 'C': 4, 'N': 2, 'O': 2}, # uracil
                    {'H': 5, 'C': 4, 'N': 3, 'O': 1}, # cytosine
                    {'H': 5, 'C': 5, 'N': 5, 'O': 1}, # guanine
                    {'H': 6, 'C': 5, 'N': 2, 'O': 2}, # thymine
                    {'H': 2, 'C': 1,         'O': 2}, # formic acid
                    {'H':10, 'C': 5,         'O': 5}] # ribose
elementsCompIndices = ['HCN', 'adenine', 'H2O', 'NH3', 'H2', 'CO', 'CH2NO',
                       'formaldehyde', 'uracil', 'cytosine', 'guanine',
                       'thymine', 'formic acid', 'ribose']
elementsComp = pd.DataFrame(elementsCompData, index=elementsCompIndices,
                            dtype=float)

# Element masses (unit [u])
elementsMass = {'H':  1.008,
                'C': 12.011,
                'N': 14.007,
                'O': 15.999}

def molar_mass(molecule):
    '''Caclulate molar mass of given molecule.

    Parameters:
        molecule : str
            Name of molecule.

    Returns:
        mass     : float
            Molar mass of molecule in unit [u] or [g/mol].
    '''
    atoms = elementsComp.loc[molecule, :]
    mass = 0
    for atom, number in atoms.items():
        if not np.isnan(number):
            mass += number * elementsMass[atom]

    return mass


def R_data(molecule, phase, pressure, tempMax):
    """Obtain data from R database 'CHNOSZ' and return it as numpy.2darray.

    Parameters:
        molecule : str
            String containing the name of the wanted molecule as stored
            in the R 'CHNOSZ' database.
        phase    : str
            String containing the name of the wanted phase of the
            molecule as stored in the R 'CHNOSZ' database.
        pressure : float
            Pressure in unit [bar].
        tempMax  : float
            Temperature up to which data should be acquired in unit
            Kelvin.

    Returns:
        out      : numpy.2darray
            Contains temperatures as first and corresponding Gibbs
            energies as second element.
    """
    tempMaxStr = '{:.2f}'.format(tempMax)

    # Use rpy2 library to start R session and obtain data
    # Load database in R
    ro.r('suppressMessages(' + 'library(CHNOSZ)' + ')')
    #ro.packages.importr('CHNOSZ', suppress_messages=True)

    # Change units to Kelvin and bar
    ro.r('suppressMessages(' + 'T.units(\'K\')' + ')')
    ro.r('suppressMessages(' + 'P.units(\'bar\')' + ')')
    ro.r('suppressMessages(' + 'E.units(\'J\')' + ')')

    # Compose command to obtain data from R as string
    Rcommand =  'subcrt(info(\'' + molecule + '\', state = \'' + phase
    Rcommand += '\')[[1]], T = seq(273.15, ' + tempMaxStr + ', '
    Rcommand += 'length.out = 200), P = ' + str(pressure) + ')'

    # Get data as R dataframe
    hyd = ro.r('suppressMessages(' + Rcommand + ')')

    # Convert to numpy array (with pandas dataframe as itermediate step)
    ro.pandas2ri.activate()
    pd_hyd = ro.conversion.rpy2py(hyd)

    out = np.transpose(pd_hyd[1][0][['T', 'G']].to_numpy())

    print('R / CHNOSZ: Obtained Gibbs energy data for \''
          + molecule
          + '\' in phase \''
          + phase
          + '\'')

    return out


def fit(molecule, phase, initConc, coeffs, pressure, tempMax, plot=True):
    """Fit Gibbs energy coefficients and return them as list of dict of list.

    Parameters:
        molecule : str
            String containing the name of the wanted molecule as stored
            in the R 'CHNOSZ' database.
        phase    : str
            String containing the name of the wanted phase of the
            molecule as stored in the R 'CHNOSZ' database.
        initConc : float
            Initial concentration in unit [mol X/mol H2O].
        coeffs   : list of dict of list
            List of dictionaries to add new fitted coefficients. First
            element contains dictionary for gas phase, second for aqeous
            phase and third for condensed phase. Keys of dictionaries are
            the names of the molecules as string and values are a list
            with numpy.array containing the coefficents as first element
            and given initial Conc 'initConc' as second element.
            Will be returned.
        pressure : float
            Pressure in unit [bar].
        tempMax  : float
            Temperature up to which data should be acquired in unit
            Kelvin.
        plot     : boolean, optional
            Whether to fit the data and the fit and store them in
            './fit_plots' subdirectory. Default is True.

    Returns:
        coeffs   : list of dict of list
            List of dictionaries to add new fitted coefficients. First
            element contains dictionary for gas phase, second for aqeous
            phase and third for condensed phase. Keys of dictionaries are
            the names of the molecules as string and values are a list
            with numpy.array containing the coefficents as first element
            and given initial concentration 'initConc' as second element.
    """
    # Fit function for Gibbs energy
    def gibbs_energy(T, a, b, c, d, e, f):
        return a + b*T + c*T*np.log(T) + d*T**2 + e*T**3 + f/T

    # Obtain data from R database 'CHNOSZ'
    data = R_data(molecule, phase, pressure, tempMax)

    # Fit gibbs energy coefficients
    popt, pcov = scipy.optimize.curve_fit(gibbs_energy, data[0], data[1])

    # Collect them in list of dict
    if phase == 'gas':
        coeffs[0][molecule] = [popt, initConc]
    elif phase == 'aq' or phase == 'liq':
        coeffs[1][molecule] = [popt, initConc]
    elif phase == 'cr':
        coeffs[2][molecule] = [popt, initConc]
    else:
        errStr  = '\'' + phase + '\' as phase argument not allowed, '
        errStr += 'has to be either \'gas\', \'aq\', \'liq\' or \'cr\'!'
        raise ValueError(errStr)

    # Plot data and fit if input parameter 'plot=True'
    if (plot):
        fig, ax = plt.subplots()

        # Data
        ax.scatter(data[0], data[1] * 1e-3, s=7,
                   label='Data from CHNOSZ database')
        ax.grid()

        # Fit
        fitLabel = 'Fit: %1.3f+%1.3fT+%1.3fTlogT+\n%1.3fT^2+%1.3fT^3+%1.3f/T' \
                   % tuple(popt)
        ax.plot(data[0], gibbs_energy(data[0], *popt) * 1e-3, color='C1',
                label=fitLabel)

        ax.set_xlabel('T [K]')
        ax.set_ylabel('Gibbs energy of formation [kJ/mol]')
        plt.title(molecule
                  + ' ('
                  + phase
                  + ') at pressure of '
                  + str(pressure)
                  + ' bar')
        plt.legend()

        # Compose file path to store plot to
        fitPlotsDir  = pathlib.Path('.') / 'fit_plots'
        # Create subdirectory if necessary
        fitPlotsDir.mkdir(exist_ok=True)
        plotFileName = (molecule
                        + '_'
                        + phase
                        + '_'
                        + str(int(pressure))
                        + 'bar_fit.pdf')
        plotPath     = fitPlotsDir / plotFileName

        # Save plot
        plt.savefig(str(plotPath))

    return coeffs


def ChemApp_file(coeffs, nucleobase, reactionNo, pressure, tempMax):
    """Write ChemApp input file and store it in './input_ChemApp/' subdirectory.
    Return initial concentrations in the order they were written to the file
    as list of float.

    Parameters:
        coeffs        : list of dict of list
            List of dictionaries containing fitted Gibbs energy
            coefficients. First element has to contain dictionary for gas
            phase, second for aqeous phase and third for condensed phase.
            Keys of dictionaries have to be the names of the molecules as
            string and values are a list with numpy.array containing the
            coefficents as first element and given initial concentration
            as second element.
        nucleobase    : str
            Name of nucleobase reaction results in.
        reactionNo    : int
            Number of reaction.
        pressure      : float
            Pressure in unit [bar].
        tempMax       : float
                 Temperature up to which data is valid in unit Kelvin.

    Returns:
        initConcs     : list of float
            List of inital concentrations of molecules in order they
            where written to ChemApp input file.
        indicesNucleo : list of int
            Index of nucleobase in initConcs.
    """
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
    fileDir.mkdir(exist_ok=True)
    filePath = fileDir / (nucleobase
                          + '_'
                          + str(reactionNo)
                          + '_'
                          + str(int(pressure))
                          + 'bar.dat')

    initConcs = []
    indicesNucleo = []
    with filePath.open(mode='w') as f:
        ####
        # Header lines
        ####
        # element1-element2-...
        for element in elementsList[:-1]:
            f.write(element + '-')
        f.write(elementsList[-1])
        today = datetime.datetime.now()
        f.write('    See commentary below!    ')
        f.write(today.strftime('MPIA %b %Y\n'))

        # Phase information
        # Number of constituents
        f.write(str(len(elementsList)) + '   ')
        # Numbers of mixing phases
        numMixPhases = 1
        for phase in coeffs[1:-1]:
            if len(phase) > 0:
                numMixPhases += 1
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
        # G(T) = a + b*T + c*T*ln(T) + d*T^2 + e*T^3 + f/T
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
                    f.write(to_FORTRAN_str('GAS') + '\n'
                            + to_FORTRAN_str('IDMX') + '\n')
                elif i == 1:
                    f.write(to_FORTRAN_str('AQUEOUS') + '\n'
                            + to_FORTRAN_str('IDMX') + '\n')

                for molecule in coeffs[i]:
                    f.write(to_FORTRAN_str(molecule.upper()) + '\n')

                    # Store initial concentrations in order they are written to
                    # file
                    initConcs.append(coeffs[i][molecule][1])

                    # Find index of nucleobase in 'initConcs'
                    if molecule == nucleobase:
                        indicesNucleo.append(int(len(initConcs) - 1))

                    f.write('1  1   ')
                    for element in elementsList:
                        f.write(str(np.nan_to_num(elementsComp.loc[molecule,
                                                                   element]
                                                 )
                                   )
                                + '   ')
                    f.write('\n')

                    # Expand string with spaces to length of 15
                    def to_15_str(string):
                        while len(string) < 15:
                            string += ' '
                        return string

                    def write_gibbs_coeffs(coeffsList, tempMax):
                        coeffsList = np.append(tempMax, coeffsList)
                        for i in range(len(coeffsList)):
                            if i%5 != 0 or i == 0:
                                valueStr = '{:.5f}'.format(coeffsList[i])
                                f.write(to_15_str(valueStr))
                            else:
                                valueStr = '{:.5f}'.format(coeffsList[i])
                                f.write('\n' + to_15_str(valueStr))

                    write_gibbs_coeffs(coeffs[i][molecule][0], tempMax)
                    f.write('\n')

        # Commentary
        f.write('\n\n########################################################')
        f.write('#######################\n\n')
        f.write('Date           : ' + today.strftime('%d.%m.%Y\n'))
        f.write('Originator     : Klaus Paschek\n')
        f.write('Serial No      :\n')
        f.write('Quality Status :\n')
        f.write('ChemSage Vers. :\n\n')
        tempMaxStr = '{:.2f}'.format(tempMax)
        f.write('Temperatures   : 273.15 to ' + tempMaxStr + ' K\n')
        f.write('Compositions   :\n')
        f.write('Applications   : Chemical equilibrium for ' + nucleobase)
        f.write(' in a meteoritic parent body.\n')
        f.write('                 Reaction No. ' + str(reactionNo) + '.\n\n\n')
        f.write('Sources of data\n\n')
        f.write('All            : CHNOSZ OrganoBioGeoTerm Database (OBIGT)\n')
        f.write('                 (http://www.chnosz.net/)\n\n')
        f.write('References     : Cobb et al. (2015),\n')
        f.write('                 Pearce and Pudritz (2016)\n\n')
        f.write('Note: This data-file models the Gibbs energy according to ')
        f.write('the equation\n\n')
        f.write('    G = a + bT + cTlnT + dT^2 + eT^3 + f/T\n\n')
        f.write('(i.e. it does not include the pressure dependence)\n\n')
        f.write('############################################################')
        f.write('###################')

    # Check if index of nucleobase was found
    if len(indicesNucleo) == 0:
        errorStr =  'Position to which nucleobase was written in ChemApp '
        errorStr += 'input file was not found'
        raise ValueError(errorStr)

    return initConcs, indicesNucleo
