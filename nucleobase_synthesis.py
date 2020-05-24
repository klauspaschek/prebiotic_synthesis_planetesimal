import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticks
import matplotlib.lines as mlines
import pathlib
import sys
import csv
import tqdm
import contextlib
import math

# Custom python libraries
import Water
import GibbsEnergy

# pybind11 C++/FORTRAN library
import ChemApp



def read_input(nucleobase, reactionNo):
    '''Read input file containing information about reactants involved in
    reaction in './reaction_info' subdirectory.

    Parameters:
        nucleobase : str
            Name of nucleobase.
        resctionNo : int
            Number of reaction.

    Returns:
        reactants  : list of list of [str, str, float]
            List containing list containing information about individual
            molecules involved in reaction. Each list contains name of molecule
            as first element, phase as second element and intial concentration
            as third element.
        waterRole  : str
            Role of water in reaction. Possible are 'solvent', 'reactant' and
            'product'.
    '''
    # Directory and file path
    inputDir      = pathlib.Path('.') / 'reaction_info'
    inputFileName = nucleobase + '_' + str(reactionNo) + '.csv'
    inputPath     = inputDir / inputFileName
    if not inputDir.is_dir():
        inputDir.mkdir()
        inputPath.open()
        errorStr =  'Directory \'./' + str(inputDir) + '\' for input file '
        errorStr += 'containing information about reactants involved in '
        errorStr += 'reaction not found'
        raise NotADirectoryError(errorStr)
    if not inputPath.is_file():
        inputPath.open()
        errorStr =  'File \'./' + str(inputPath) + '\' '
        errorStr += 'containing information about reactants involved in '
        errorStr += 'reaction not found'
        raise FileNotFoundError(errorStr)

    # Read
    reactants = []
    with inputPath.open() as f:
        readReactants = csv.reader(f, delimiter=',')
        for row in readReactants:
            reactants.append(row)

    # Find role of water in reaction
    waterRole = reactants[0][1]
    if waterRole != 'solvent' and waterRole != 'reactant' \
       and waterRole != 'product':
        errorStr =  'Role of water in reaction has to be specified either as '
        errorStr += '\'solvent\', \'reactant\' or \'product\' as second '
        errorStr += 'element in first row of file \''
        errorStr += str(inputPath)
        erorrStr += '\'. E.g.:\nH2O,solvent'
        raise ValueError(errorStr)

    # Select information about reactants
    reactants = reactants[1:]

    # Convert last element of each list containing intial concentration in
    # 'reactants' to float
    for i in range(len(reactants)):
        reactants[i][2] = float(reactants[i][2])

    return reactants, waterRole

def write_input_ChemApp(nucleobase, reactionNo, pressure, reactants):
    '''Create thermochemical input file for ChemApp.

    Parameters:
        nucleobase    : str
            Name of nucleobase.
        reactionNo    : int
            Number of reaction.
        pressure      : float
            Pressure in unit [bar].
        reactants     : list of list of [str, str, float]
            List containing list containing information about individual
            molecules involved in reaction. Each list contains name of molecule
            as first element, phase as second element and intial concentration
            as third element.

    Returns:
        initConcs     : list of float
            List of inital concentrations of molecules in order they
            where written to ChemApp input file.
        indicesNucleo : list of int
            Index of nucleobase in initConcs.
        tempMax       : float
            Temperature up to which data is valid in unit Kelvin.
    '''
    # Calculate temperature of boiling point of water as upper bound for data
    # acquisition in R library 'CHNOSZ'
    tempMax = Water.boiling_point(pressure)
    # Calculate Gibbs energy coefficients for reactants read from input file
    # using data from R library 'CHNOSZ'
    coeffs = [{}, {}, {}]
    for reactant in reactants:
        coeffs = GibbsEnergy.fit(reactant[0], reactant[1], reactant[2], coeffs,
                                 pressure, tempMax)

    # Write input file for ChemApp and obtain initial concentration in correct
    # order
    initConcs, indicesNucleo = GibbsEnergy.ChemApp_file(coeffs,
                                                        nucleobase,
                                                        reactionNo,
                                                        pressure, tempMax)

    return initConcs, indicesNucleo, tempMax

def iter_temps_amounts(nucleobase, reactionNo, pressure, temps,
                       recursion=False, numbering=False):
    '''Iterate over given temperatures in parameter 'temps' and return
    corresponding equilibrium amounts of nucleobase.

    Parameters:
        nucleobase : str
            Name of nucleobase.
        reactionNo : int
            Number of reaction.
        pressure   : float
            Pressure in unit [bar].
        temps      : numpy.array of float
            Temperatures to calculate equilibrium amounts for in unit [Kelvin].
        recursion  : boolean, optional
            If True, amounts of constituents of previous iteration are used
            for next iteration as initial concentrations.
        numbering  : boolena, optional
            If True, write a separate output '.csv' file in './results'
            subdirectory for each call of this function with increasing number
            in file name.

    Returns:
        amounts    : numpy.array of float
            Calculated equilibrium amounts corresponding to given temperatures
            in 'temps' in unit [mol nucleobase/mol water].
    '''
    # Obtain information about reaction from input file
    reactants, waterRole = read_input(nucleobase, reactionNo)

    # Write input file for ChemApp
    initConcs, indicesNucleo, tempMax = write_input_ChemApp(nucleobase,
                                                            reactionNo,
                                                            pressure,
                                                            reactants)

    # Error variable to check execution of FORTRAN routines in ChemApp
    err = int(0)

    # Directory for ChemApp to write log files to
    logDir = pathlib.Path('.') / 'logs'
    logDir.mkdir(exist_ok=True)

    # Path of thermochemical input file for ChemApp
    fileDir = pathlib.Path('.') / 'input_ChemApp'
    if not fileDir.is_dir():
        errorStr =  'Directory \'' + str(fileDir) + '\'necessary for ChemApp '
        errorStr += 'to read input file from not found'
        raise NotADirectoryError(errorStr)
    filePath = fileDir / (nucleobase
                          + '_'
                          + str(reactionNo)
                          + '_'
                          + str(int(pressure))
                          + 'bar.dat')
    if not filePath.is_file():
        errorStr =  'File \'./' + str(filePath) + '\' '
        errorStr += 'containing thermochemical information about reactants '
        errorStr += 'involved in reaction for ChemApp to read from not found'
        raise FileNotFoundError(errorStr)

    # Read thermochemical input file into ChemApp
    ChemApp.read_data(str(filePath), err)

    # Do iteration over different temperatures
    amounts = np.zeros(len(temps))
    result = initConcs.copy()
    for i in range(len(temps)):
        # Check if water is liquid at this temperature to allow synthesis of
        # nucleobase
        if not recursion:
            if temps[i] >= 273.15 and temps[i] <= tempMax: # liquid water
                result = ChemApp.amounts(initConcs, waterRole,
                                         temps[i], pressure,
                                         str(logDir), err)
                # Sum up all constituents of nucleobase
                amountNucleo = 0
                for index in indicesNucleo:
                    amountNucleo += result[index]
                amounts[i] = amountNucleo
            else:
                amounts[i] = 0.0#np.nan

        else: # if recursion
            if temps[i] >= 273.15 and temps[i] <= tempMax: # liquid water
                result = ChemApp.amounts(result, waterRole,
                                         temps[i], pressure,
                                         str(logDir), err)
                # Sum up all constituents of nucleobase
                amountNucleo = 0
                for index in indicesNucleo:
                    amountNucleo += result[index]
                amounts[i] = amountNucleo
            elif temps[i] <= tempMax: # frozen, formed nucleobase is preserved
                amounts[i] = amounts[i-1]
            else: # too hot, water boils and destroys formed nucleobase
                amounts[i] = 0.0
                for index in indicesNucleo:
                    result[index] = 0.0#amounts[i]

    # Save amounts in .csv file
    # Compose file path to store .csv file to
    amountsDir = pathlib.Path('.') / 'results'
    # Create subdirectory if necessary
    amountsDir.mkdir(exist_ok=True)
    if numbering:
        amountsFile = amountsDir / (nucleobase
                                    + '_'
                                    + str(reactionNo)
                                    + '_'
                                    + str(int(pressure))
                                    + 'bar_amounts_'
                                    + str(iter_temps_amounts.counter)
                                    + '.csv')
        iter_temps_amounts.counter += 1
    else:
        amountsFile = amountsDir / (nucleobase
                                    + '_'
                                    + str(reactionNo)
                                    + '_'
                                    + str(int(pressure))
                                    + 'bar_amounts.csv')

    # Write .csv file
    with amountsFile.open(mode='w') as f:
        # Write header
        f.write('Temperature[K],Amount[mol ' + nucleobase + '/mol H2O]\n')

        for i in range(len(temps)):
            f.write('{:.2f},{:.5e}\n'.format(temps[i], amounts[i]))

    return amounts

def sci_fmt(x, pos):
    if x <= 0:
        return u"${:.0f}$".format(x)
    else:
        f = mticks.ScalarFormatter(useOffset=False, useMathText=True)
        return u"${}$".format(f._formatSciNotation('%1.10e' % x))

def plot_peak_temps(tempsData, targetNucleobases, reactionNos, pressure,
                    rhoI, rhoP, phi, tOF):
    radii = tempsData[1:, 0] * 1e-3 # unit [km]
    temps = np.amax(tempsData[1:, 1:], axis = 1) # unit [Kelvin]

    # Error variable to check execution of FORTRAN routines in ChemApp
    err = int(0)

    # Startup ChemApp
    ChemApp.start(err)

    # Directory for ChemApp to write log files to
    logDir = pathlib.Path('.') / 'logs'
    logDir.mkdir(exist_ok=True)

    # Open log files
    ChemApp.open_file(str(logDir), err);

    # Calculate amounts of nucleobase(s)
    amounts = []
    with std_out_err_redirect_tqdm() as origStdOut:
        with tqdm.tqdm(total=sum([len(nos) for nos in reactionNos]),
                       file=origStdOut, dynamic_ncols=True) as pbar:
            for i in range(len(targetNucleobases)):
                amountsNucleobase = []
                for j in range(len(reactionNos[i])):
                    amountsNucleobaseNo = iter_temps_amounts(
                                              targetNucleobases[i],
                                              reactionNos[i][j],
                                              pressure,
                                              temps)

                    # Unit conversion
                    # from [mol nucleobase/mol water]
                    # to [kg nucleobase/kg planetesimal] in [ppb]
                    amountsNucleobaseNo *= ((GibbsEnergy.molar_mass(
                                                 targetNucleobases[i])
                                             / GibbsEnergy.molar_mass('H2O'))
                                            * (rhoI / rhoP)
                                            * phi
                                            * 1e9)

                    amountsNucleobase.append(amountsNucleobaseNo)
                    pbar.update(1)
                amounts.append(amountsNucleobase)

    # Close log files
    ChemApp.close_file(err)

    # Plotting
    # Styling
    fontSize = 8
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    del colors[7]
    del colors[5]
    del colors[4]
    del colors[1]
    linestyles = [(0, (2, 12)), (2, (2, 12)), (4, (2, 12)), (6, (2, 12)),
                  (8, (2, 12)), (10, (2, 12))]
    markers = [None, 'x']
    markeverys = [None, (0, 9), None, (3, 9), None, (6, 9)]
    legendHandles=[]
    styleIndex = 0

    # Plot amounts
    fig, ax1 = plt.subplots()
    for i in range(len(targetNucleobases)):
        for j in range(len(reactionNos[i])):
            ax1.plot(radii, amounts[i][j],
                     color=colors[styleIndex % len(colors)],
                     linestyle=linestyles[styleIndex % len(linestyles)],
                     marker=markers[styleIndex % len(markers)],
                     markevery=markeverys[styleIndex % len(markeverys)])
            legendHandles.append(mlines.Line2D([], [],
                                 color=colors[styleIndex % len(colors)],
                                 marker=markers[styleIndex % len(markers)],
                                 label=plotLabelIndex[reactionNos[i][j]]))
            styleIndex += 1

    ax1.set_xlim(np.min(radii), math.ceil(np.max(radii)))
    ax1.set_xlabel('Radius [km]', fontsize=fontSize)
    ax1.set_yscale('log')
    ax1.set_ylim(1e0, 3e6)
    ax1.set_ylabel('Nucleobase abundance [ppb]', fontsize=fontSize)
    ax1.tick_params(axis='both', which='both', top=True, direction='in',
                    labelsize=fontSize)

    ax1.legend(handles=legendHandles, prop={'size': fontSize})

    # Plot temperatures
    color = 'C1'
    ax2 = ax1.twinx()
    ax2.plot(radii, temps, color=color)
    ax2.set_ylabel(r'T$_{max}$ [K]', color=color, fontsize=fontSize)
    ax2.tick_params(axis='y', direction='in', labelcolor=color,
                    labelsize=fontSize)

    plt.title('Distribution of nucleobase abundances \ndepending on distance '
              + 'from center of planetesimal at peak temperature.\nProperties of planetesimal: '
              + r'Radius = '
              + str(math.ceil(np.max(radii)))
              + r' km, $\rho$ = '
              + str(rhoP)
              + r' kg/m$^3$, $\rho_{ice}$ = '
              + str(rhoI)
              + r' kg/m$^3$, porosity = '
              + str(phi)
              + ',\ntime of formation after CAI = '
              + str(tOF)
              + ' Myr.',
              fontdict = {'fontsize': fontSize})

    nucleobaseTextStr = ''
    if len(targetNucleobases) == 1:
        nucleobaseTextStr = (targetNucleobases[0][0].upper()
                             + ' Synthesis')
    elif len(targetNucleobases) == 2:
        nucleobaseTextStr = (targetNucleobases[0][0].upper()
                             + ' and '
                             + targetNucleobases[1][0].upper()
                             + ' Synthesis')
    else:
        for nucleobase in targetNucleobases[:-1]:
            nucleobaseTextStr += str(nucleobase[0].upper()) + ', '
        nucleobaseTextStr = nucleobaseTextStr[:-2]
        nucleobaseTextStr += (' and '
                              + targetNucleobases[-1][0].upper()
                              + ' Synthesis')

    ax2.text(0.02, 0.96, nucleobaseTextStr, transform=ax1.transAxes,
             fontsize=fontSize, fontweight='bold')

    fig.tight_layout()

    # Save plot
    plotDir = pathlib.Path('.') / 'results'
    plotDir.mkdir(exist_ok=True)
    nucleobaseFileNameStr = ''
    for nucleobase in targetNucleobases:
        nucleobaseFileNameStr += nucleobase + '_'
    plotPath = plotDir / (nucleobaseFileNameStr
                          + str(int(pressure))
                          + 'bar_peak_temps_amounts.pdf')

    plt.savefig(str(plotPath), bbox_inches = 'tight', pad_inches = 0.01)

def plot_time_iter(tempsData, targetNucleobases, reactionNo, pressure, step,
                   rhoI, rhoP, phi, tOF):
    time  = tempsData[0, 1::step] * 1e-6  # unit [Myr]
    temps = tempsData[1, 1::step]         # unit [K]

    # Error variable to check execution of FORTRAN routines in ChemApp
    err = int(0)

    # Startup ChemApp
    ChemApp.start(err)

    # Directory for ChemApp to write log files to
    logDir = pathlib.Path('.') / 'logs'
    logDir.mkdir(exist_ok=True)

    # Open log files
    ChemApp.open_file(str(logDir), err);

    amounts = []
    iter_temps_amounts.counter = 0
    with std_out_err_redirect_tqdm() as origStdOut:
        with tqdm.tqdm(total=sum([len(nos) for nos in reactionNos]),
                       file=origStdOut, dynamic_ncols=True) as pbar:
            for i in range(len(targetNucleobases)):
                amountsNucleobase = []
                for j in range(len(reactionNos[i])):
                    amountsNucleobaseNo = iter_temps_amounts(
                                              targetNucleobases[i],
                                              reactionNos[i][j],
                                              pressure,
                                              temps,
                                              recursion=True, numbering=True)

                    # Unit conversion
                    # from [mol nucleobase/mol water]
                    # to [kg nucleobase/kg planetesimal] in [ppb]
                    amountsNucleobaseNo *= ((GibbsEnergy.molar_mass(
                                                 targetNucleobases[i])
                                             / GibbsEnergy.molar_mass('H2O'))
                                            * (rhoI / rhoP)
                                            * phi
                                            * 1e9)

                    amountsNucleobase.append(amountsNucleobaseNo)
                    pbar.update(1)
                amounts.append(amountsNucleobase)

    # Close log files
    ChemApp.close_file(err)

    # Plotting
    # Styling
    fontSize = 8
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    del colors[7]
    del colors[5]
    del colors[4]
    del colors[1]
    linestyles = [(0, (2, 12)), (2, (2, 12)), (4, (2, 12)), (6, (2, 12)),
                  (8, (2, 12)), (10, (2, 12))]
    markers = [None, 'x']
    markeverys = [None, (0, 3), None, (1, 3), None, (1, 3)]
    legendHandles=[]
    styleIndex = 0

    # Plot amounts
    fig, ax1 = plt.subplots()
    for i in range(len(targetNucleobases)):
        for j in range(len(reactionNos[i])):
            ax1.plot(time, amounts[i][j],
                     color=colors[styleIndex % len(colors)],
                     linestyle=linestyles[styleIndex % len(linestyles)],
                     marker=markers[styleIndex % len(markers)],
                     markevery=markeverys[styleIndex % len(markeverys)])
            legendHandles.append(mlines.Line2D([], [],
                                 color=colors[styleIndex % len(colors)],
                                 marker=markers[styleIndex % len(markers)],
                                 label=plotLabelIndex[reactionNos[i][j]]))
            styleIndex += 1

    ax1.set_xscale('log')
    ax1.set_xlabel('Time after formation [Myr]', fontsize=fontSize)
    ax1.set_xlim(1e-3, 5e3)
    ax1.set_yscale('log')
    ax1.set_ylim(1e0, 3e6)
    ax1.set_ylabel('Nucleobase abundance [ppb]', fontsize=fontSize)
    ax1.tick_params(axis='both', which='both', top=True, direction='in',
                    labelsize=fontSize)

    # Plot temperatures
    color = 'C1'
    ax2 = ax1.twinx()
    ax2.plot(time, temps, color=color)
    ax2.set_ylabel(r'T$_{core}$ [K]', color=color, fontsize=fontSize)
    ax2.tick_params(axis='y', direction='in', labelcolor=color,
                    labelsize=fontSize)

    radii = tempsData[1:, 0] * 1e-3
    plt.title('Temporal evolution of nucleobase abundances in the center of '
              + 'the planetesimal.\nProperties of planetesimal: '
              + r'Radius = '
              + str(math.ceil(np.max(radii)))
              + r' km, $\rho$ = '
              + str(rhoP)
              + r' kg/m$^3$, $\rho_{ice}$ = '
              + str(rhoI)
              + r' kg/m$^3$, porosity = '
              + str(phi)
              + ',\ntime of formation after CAI = '
              + str(tOF)
              + ' Myr.',
              fontdict = {'fontsize': fontSize})

    nucleobaseTextStr = ''
    if len(targetNucleobases) == 1:
        nucleobaseTextStr = (targetNucleobases[0][0].upper()
                             + ' Synthesis')
    elif len(targetNucleobases) == 2:
        nucleobaseTextStr = (targetNucleobases[0][0].upper()
                             + ' and '
                             + targetNucleobases[1][0].upper()
                             + ' Synthesis')
    else:
        for nucleobase in targetNucleobases[:-1]:
            nucleobaseTextStr += str(nucleobase[0].upper()) + ', '
        nucleobaseTextStr = nucleobaseTextStr[:-2]
        nucleobaseTextStr += (' and '
                              + targetNucleobases[-1][0].upper()
                              + ' Synthesis')

    ax1.text(0.02, 0.96, nucleobaseTextStr, transform=ax1.transAxes,
             fontsize=fontSize, fontweight='bold')

    ax1.legend(handles=legendHandles, prop={'size': fontSize},
               loc='best', bbox_to_anchor=(0., 0.55, 1., 0.41))

    fig.tight_layout()

    # Save plot
    plotDir = pathlib.Path('.') / 'results'
    plotDir.mkdir(exist_ok=True)
    nucleobaseFileNameStr = ''
    for nucleobase in targetNucleobases:
        nucleobaseFileNameStr += nucleobase + '_'
    plotPath = plotDir / (nucleobaseFileNameStr
                          + str(int(pressure))
                          + 'bar_time_iter_amounts.pdf')

    plt.savefig(str(plotPath), bbox_inches = 'tight', pad_inches = 0.01)

'''
def plot_time_iter_radii(tempsData, targetNucleobase, reactionNo, pressure,
                         step,
                         rhoI, rhoP, phi):
    radii = tempsData[1::10, 0      ]         # unit [km]
    time  = tempsData[0    , 1::step] * 1e-6  # unit [Myr]
    temps = tempsData[1::10, 1::step]         # unit [K]

    # Error variable to check execution of FORTRAN routines in ChemApp
    err = int(0)

    # Startup ChemApp
    ChemApp.start(err)

    # Directory for ChemApp to write log files to
    logDir = pathlib.Path('.') / 'logs'
    logDir.mkdir(exist_ok=True)

    # Open log files
    ChemApp.open_file(str(logDir), err);

    amounts = np.zeros([len(temps), len(time)])
    iter_temps_amounts.counter = 0
    with std_out_err_redirect_tqdm() as origStdOut:
        for i in tqdm.tqdm(range(len(temps)), file=origStdOut,
                           dynamic_ncols=True):
            amounts[i] = iter_temps_amounts(targetNucleobase, reactionNo,
                                            pressure,
                                            temps[i],
                                            recursion=True, numbering=True)

    # Close log files
    ChemApp.close_file(err)

    # Unit conversion from
    # [mol nucleobase/mol water] to [kg nucleobase/kg planetesimal] in [ppb]
    amounts *= ((GibbsEnergy.molar_mass(targetNucleobase)
                 / GibbsEnergy.molar_mass('H2O'))
                * (rhoI / rhoP)
                * phi
                * 1e9)

    # Plot amounts
    fig, ax = plt.subplots()
    for i in range(len(amounts)):
        ax.plot(time, amounts[i], label='{:2.0f} km'.format(radii[i]))
    ax.set_xscale('log')
    ax.set_xlabel('Time after formation [Myr]')
    ax.set_ylabel('Nucleobase abundance [ppb]')
    ax.set_xlim(1e-3, 5e3)

    #ax.yaxis.set_major_formatter(mtick.FuncFormatter(sci_fmt))

    plt.title('Temporal evolution of nucleobase $\\bf{'
              + targetNucleobase
              + '}$ in reaction no. $\\bf{'
              + str(reactionNo)
              + '}$\ndepending on distance from center of planetesimal',
              fontdict = {'fontsize': 10})
    plt.legend()
    fig.tight_layout()

    # Save plot
    plotDir = pathlib.Path('.') / 'results'
    plotDir.mkdir(exist_ok=True)
    plotPath = plotDir / (targetNucleobase
                          + '_'
                          + str(reactionNo)
                          + '_'
                          + str(int(pressure))
                          + 'bar_time_iter_amounts.pdf')

    plt.savefig(str(plotPath))
'''



####
# Redirect print to 'tqdm.tqdm.write()' for correct behavior of progress bar
####

class DummyTqdmFile(object):
    """Dummy file-like that will write to tqdm"""
    file = None
    def __init__(self, file):
        self.file = file
    def flush(self):
        getattr(self.file, 'flush', lambda: None)()
    def write(self, x):
        # Avoid print() second call (useless \n)
        if len(x.rstrip()) > 0:
            tqdm.tqdm.write(x, file=self.file)

@contextlib.contextmanager
def std_out_err_redirect_tqdm():
    origOutErr = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = map(DummyTqdmFile, origOutErr)
        yield origOutErr[0]
    # Relay exceptions
    except Exception as exc:
        raise exc
    # Always restore sys.stdout/err if necessary
    finally:
        sys.stdout, sys.stderr = origOutErr



####
# --------------------------- Main Program ------------------------------------
####

####
# Command line arguments
####

# Dictionary with all possible reaction numbers
reactionIndex = {'adenine' : [1, 3, 4, 6, 7, 8],
                 'uracil'  : [29, 32],
                 'cytosine': [43, 44],
                 'guanine' : [51, 54],
                 'thymine' : [58, 62]}

plotLabelIndex = {1: r'1. CO + H$_2$ + NH$_3$ $\longrightarrow$ Adenine '
                     + '+ H$_2$O',
                  3: r'3. 5HCN$_{(aq)}$ $\longrightarrow$ Adenine$_{(aq)}$',
                  4: r'4. HCN + NH$_3$ $\longrightarrow$ Adenine',
                  6: r'6. 5CO + 5NH$_3$ $\longrightarrow$ Adenine '
                     + '+ 5H$_2$O',
                  7: r'7. HCN + H$_2$O $\longrightarrow$ Adenine',
                  8: r'8. HCN + NH$_3$ + H$_2$O $\longrightarrow$ Adenine',
                  29: r'29. 2HCN$_{(aq)}$ + 2CH$_2$O$_{(aq)}$ '
                      + r'$\longrightarrow$ Uracil$_{(aq)}$ + H$_{2(aq)}$',
                  32: r'32. Cytosine + H$_2$O $\longrightarrow$ Uracil + '
                      + r'NH$_3$',
                  43: r'43. CO + H$_2$ + NH$_3$ $\longrightarrow$ Cytosine + '
                      + r'H$_2$O',
                  44: r'44. 3HCN$_{(aq)}$ + CH$_2$O$_{(aq)}$ $\longrightarrow$'
                      + r'Cytosine$_{(aq)}$',
                  51: r'51. CO + H$_2$ + NH$_3$ $\longrightarrow$ Guanine + '
                      + r'H$_2$O',
                  54: r'54. 5HCN$_{(aq)}$ + H$_2$O $\longrightarrow$ '
                      + r'Guanine$_{(aq)}$ + H$_2(aq)$',
                  58: r'58. 2HCN$_{(aq)}$ + 3CH$_2$O$_{(aq)}$ '
                      + r'$\longrightarrow$ Thymine$_{(aq)}$ + H$_2$O',
                  62: r'62. Uracil + CH$_2$O + CH$_2$O$_2$ + H$_2$O '
                      + r'$\longrightarrow$ Thymine'}

# Check if command line arguments are provided
if len(sys.argv) < 2:
    errorStr =  'Command line argument(s) missing. When running '
    errorStr += 'this program, you have to provide at least 1 command line '
    errorStr += 'argument:\npython3 nucleobase_synthesis.py <target '
    errorStr += 'nucleobase>[/target nucleobase2>/...] [<reaction no.> or all]'
    errorStr += ' [<pressure in bar>]'
    raise ValueError(errorStr)

# Assign command line arguments for later use
targetNucleobases = sys.argv[1].split('/')
if len(targetNucleobases) == 1:
    if len(sys.argv) == 2:
        reactionNos = [reactionIndex[targetNucleobases[0]]]
        pressure = float(100)
    elif len(sys.argv) == 3:
        if sys.argv[2] == 'all':
            reactionNos = [reactionIndex[targetNucleobases[0]]]
        else:
            reactionNos = [[int(sys.argv[2])]]
        pressure = float(100)
    elif len(sys.argv) == 4:
        if sys.argv[2] == 'all':
            reactionNos = [reactionIndex[targetNucleobases[0]]]
        else:
            reactionNos = [[int(sys.argv[2])]]
        pressure = float(sys.argv[3])
    elif len(sys.argv) > 4:
        errorStr = 'Too many command line arguments provided'
        raise ValueError(errorStr)
else:
    reactionNos = []
    for nucleobase in targetNucleobases:
        reactionNos.append(reactionIndex[nucleobase])
    if len(sys.argv) == 2:
        pressure = float(100)
    elif len(sys.argv) == 3:
        pressure = float(sys.argv[2])
    else:
        errorStr = 'Too many command line arguments provided'
        raise ValueError(errorStr)


####
# Read temperature input file
####

tempsDir  = pathlib.Path('.') / 'temps_input'
tempsPath = tempsDir          / '100km.csv'

tempsData = np.genfromtxt(tempsPath, delimiter=',', unpack=True)

####
# Define properties of asteroid
####

rhoI = 917  # density of water ice [kg/m^3]
rhoP = 3000 # density of planetesimal [kg/m^3]
phi  = 0.2  # porosity
tOF  = 3.5  # time of formation after CAI [Myr]

####
# Plot
####

# Disable warning of matplotlib for too many open figures
plt.rcParams.update({'figure.max_open_warning': 0})

plot_peak_temps(tempsData, targetNucleobases, reactionNos, pressure,
                rhoI, rhoP, phi, tOF)

step = 1
plot_time_iter(tempsData, targetNucleobases, reactionNos, pressure, step,
               rhoI, rhoP, phi, tOF)

print('\n*** DONE ***')
