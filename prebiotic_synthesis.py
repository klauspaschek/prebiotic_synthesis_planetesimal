import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticks
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
# Use desired matplotlib stylefile for setting fonts
plt.style.use('./latex.mplstyle')
# Load LaTeX libraries for matplotlib
params = {'text.latex.preamble' : [r'\usepackage{mhchem}']}#, palatino, mathpazo}']}
plt.rcParams.update(params)
import pathlib
import sys
import csv
import tqdm
import contextlib
import math
import warnings

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
                       recursion=False, numbering=False, waterConc=None):
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

    if not waterConc == None:
        for i in range(len(initConcs)):
            if not initConcs[i] == 1.0:
                initConcs[i] = initConcs[i] * waterConc
    #initConcs = [i * waterConc for i in initConcs]

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


####
# Helper functions for example plotting routines below
####

def sci_fmt(x, pos):
    if x <= 0:
        return u"${:.0f}$".format(x)
    else:
        f = mticks.ScalarFormatter(useOffset=False, useMathText=True)
        return u"${}$".format(f._formatSciNotation('%1.10e' % x))

def token(nucleobase):
    if nucleobase == 'C128':
        return nucleobase
    elif nucleobase == 'adenosine':
        return 'A-R'
    elif nucleobase == 'uridine':
        return 'U-R'
    elif nucleobase == 'cytidine':
        return 'C-R'
    elif nucleobase == 'guanosine':
        return 'G-R'
    elif nucleobase == 'thymidine':
        return 'T-R'
    else:
        return nucleobase[0].upper()

def compValues(targetNubleobases):
    result = []
    for nucleobase in targetNucleobases:
        if nucleobase == 'adenine':
            result.append([1, 267, 'CM2 A', 'salmon'])
            #result.append([1, 267, 'CM2 Adenine Abundances', 'salmon'])
        elif nucleobase == 'guanine':
            result.append([21, 515, 'CM2 G', 'goldenrod'])
            #result.append([21, 515, 'CM2 Guanine Abundances', 'goldenrod'])
        elif nucleobase == 'uracil':
            result.append([37, 63, 'CM2 U', 'lawngreen'])
            #result.append([37, 63, 'CM2 Uracil Abundances', 'lawngreen'])
        elif nucleobase == 'ribose':
            result.append([4.5, 25, 'CM2/CR2 Ribose Abundances',
                           'lightskyblue'])

    return result



####
# Example plotting routines
####

def plot_peak_temps(tempsData, targetNucleobases, reactionNos, pressure,
                    rhoI, rhoP, phi, tOF, origStdOut, waterConc=None):
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
    with tqdm.tqdm(total=sum([len(nos) for nos in reactionNos]),
                   desc=' 1 of 2',
                   leave=False, position=0,
                   file=origStdOut, dynamic_ncols=True) as pbar:
        for i in range(len(targetNucleobases)):
            amountsNucleobase = []
            for j in range(len(reactionNos[i])):
                amountsNucleobaseNo = iter_temps_amounts(
                                          targetNucleobases[i],
                                          reactionNos[i][j],
                                          pressure,
                                          temps,
                                          waterConc=waterConc)

#                if targetNucleobases[i] == 'ribose':
#                    amountsNucleobaseNo *= yieldRibose

                # Unit conversion
                # from [mol nucleobase/mol water]
                # to [kg nucleobase/kg planetesimal] in [ppb]
                amountsNucleobaseNo *= ((GibbsEnergy.molar_mass(
                                             targetNucleobases[i])
                                         / GibbsEnergy.molar_mass('H2O'))
                                        * (rhoI / rhoP)
                                        * phi
                                        * 1e9)

                if not targetNucleobases[i] == 'ribose':
                    amountsNucleobase.append(amountsNucleobaseNo)
                else:
                    for catalyst in range(len(yieldRibose)):
                        amountsNucleobase.append(amountsNucleobaseNo
                                                 * yieldRibose[catalyst][1])

                pbar.update(1)
            amounts.append(amountsNucleobase)

    # Close log files
    ChemApp.close_file(err)

    # Plotting
    # Styling
    #fontSize = 10
    plt.rcParams["figure.figsize"] = (7.101*0.5, 7.101*0.5/1.61803)
    linestyles = [(0, (2, 10)), (2, (2, 10)), (4, (2, 10)), (6, (2, 10)),
                  (8, (2, 10)), (10, (2, 10))]
    #linestyles = [(0, (1, 6)),
    #              (1, (1, 5)),
    #              (2, (1, 4)),
    #              (3, (1, 3)),
    #              (4, (1, 2)),
    #              (5, (1, 1))]
    colors = plt.cm.viridis(np.linspace(1/6, 5/6, 5))
    colors = np.concatenate(([[0.01, 0.01, 0.01, 1.]], colors))
    markers = ['x', 'd', 'o', '^', 's', 'p']
    markeverys = [(0, 18), (3, 18), (6, 18), (9, 18), (12, 18), (15, 18)]
    legendHandles=[]
    styleIndex = 0


    # Plot amounts
    fig, ax1 = plt.subplots()
    for i in range(len(targetNucleobases)):
        if not targetNucleobases[i] == 'ribose':
            for j in range(len(reactionNos[i])):
                ax1.plot(radii, amounts[i][j],
                         color=colors[styleIndex % len(colors)],
                         linestyle=linestyles[styleIndex % len(linestyles)],
                         marker=markers[styleIndex % len(markers)],
                         markevery=markeverys[styleIndex % len(markeverys)],
                         markersize=5)
                legendHandles.append(mlines.Line2D([], [],
                                     color=colors[styleIndex % len(colors)],
                                     marker=markers[styleIndex % len(markers)],
                                     markersize=5,
                                     label=plotLabelIndex[reactionNos[i][j]]))
                styleIndex += 1
        else:
            for j in range(len(yieldRibose)):
                ax1.plot(radii, amounts[i][j],
                         color=colors[styleIndex % len(colors)],
                         linestyle=linestyles[styleIndex % len(linestyles)],
                         marker=markers[styleIndex % len(markers)],
                         markevery=markeverys[styleIndex % len(markeverys)],
                         markersize=5)
                legendHandles.append(mlines.Line2D([], [],
                                     color=colors[styleIndex % len(colors)],
                                     marker=markers[styleIndex % len(markers)],
                                     markersize=5,
                                     label=plotLabelIndex[yieldRibose[j][0]]))
                styleIndex += 1

    # Shade regions representing experimentally found abundances in
    # C chondrites for comparision
    comp = compValues(targetNucleobases)
    for elem in comp:
        ax1.axhspan(elem[0], elem[1], alpha=0.3, color=plt.cm.viridis(0.99))#elem[3])
        legendHandles.append(mpatches.Patch(label=elem[2],
                                            color=plt.cm.viridis(0.99),#elem[3],
                                            alpha=0.3))

    ax1.set_xlim(np.min(radii), math.ceil(np.max(radii)))
    ax1.set_xlabel('Radius [km]')#, fontsize=fontSize)
    ax1.set_yscale('log')
    #_, ylimTop = ax1.get_ylim()
    ax1.set_ylim(1e0, 6e2)
    ax1.set_ylabel('Molecular abundance [ppb]')#, fontsize=fontSize)
    ax1.tick_params(axis='both', which='both', top=True, direction='in')#,
#                    labelsize=fontSize)

    # Plot temperatures
    color = plt.cm.viridis(0)
    ax2 = ax1.twinx()
    ax2.plot(radii, temps, color=color)
    legendHandles.append(mlines.Line2D([], [],
                         color=color,
                         label=r'Temperature $T_{\mathrm{max}}$'))
    ax2.set_xlim(0,150)
    #_, ylimTop = ax2.get_ylim()
    ax2.set_ylim(200, 525)#ylimTop * 1.03)
    ax2.set_ylabel(r'$T_{\mathrm{max}}$ [K]', color=color)#, fontsize=fontSize)
    ax2.tick_params(axis='y', direction='in', labelcolor=color)#,
#                    labelsize=fontSize)
    ax2.set_zorder(ax1.get_zorder() + 1)
    ax1.patch.set_visible(False)

#    titleStr = ('Distribution of molecular abundances \ndepending on distance '
#                + 'from center of planetesimal at peak temperature.\n'
#    titleStr = ('Properties of planetesimal: Radius = '
#                + '{}'.format(int(math.ceil(np.max(radii))))
#                + r' km, $\rho_{\mathrm{rock}}$ = '
#                + str(rhoP)
#                + r' kg/m$^3$, $\rho_{\mathrm{ice}}$ = '
#                + str(rhoI)
#                + r' kg/m$^3$, porosity = '
#                + str(phi)
#                + ',\ntime of formation after CAI = '
#                + str(tOF)
#                + ' Myr.')
#    if 'ribose' in targetNucleobases:
#        titleStr += ('\nExperimentally found yield of ribose within 5C sugars '
#                     + 'used: '
#                     + str(yieldRibose)
#                     + '.')
#    if not waterConc == None:
#        titleStr += ('\nInitial water concentration changed by factor of '
#                     + str(waterConc)
#                     + '.')
#    plt.title(titleStr, fontdict = {'fontsize': fontSize})

    # Create text label to give overview over targeted molecules
    nucleobaseTextStr = ''
    if len(targetNucleobases) == 1:
        nucleobaseTextStr = (token(targetNucleobases[0])
                             + ' Synthesis')
    elif len(targetNucleobases) == 2:
        nucleobaseTextStr = (token(targetNucleobases[0])
                             + ' and '
                             + token(targetNucleobases[1])
                             + ' Synthesis')
    else:
        for nucleobase in targetNucleobases[:-1]:
            nucleobaseTextStr += (token(nucleobase) + ', ')
        nucleobaseTextStr = nucleobaseTextStr[:-2]
        nucleobaseTextStr += (' and '
                              + token(targetNucleobases[-1])
                              + ' Synthesis')

    ax1.text(0.025, 0.915, r'\textbf{' + nucleobaseTextStr + '}',
             transform=ax1.transAxes)#, fontsize=fontSize)#, fontweight='bold'

    # Legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        ax2.legend(lines1 + lines2, labels1 + labels2, handles=legendHandles,
                   loc=4, prop={'size': 5},#-3},
                   bbox_to_anchor=(0.01, 0.01, 0.98, 0.86))

    fig.tight_layout(pad=0)

    # Save plot
    plotDir = pathlib.Path('.') / 'results'
    plotDir.mkdir(exist_ok=True)
    nucleobaseFileNameStr = ''
    for nucleobase in targetNucleobases:
        nucleobaseFileNameStr += nucleobase + '_'
    plotPath = plotDir / (nucleobaseFileNameStr
                          + str(int(pressure))
                          + 'bar_peak_temps_amounts_radius_'
                          + '{}'.format(int(math.ceil(np.max(radii))))
                          + 'km.pdf')

    plt.savefig(str(plotPath),# bbox_inches = 'tight', pad_inches = 0,
                transparent=True)


def plot_time_iter(tempsData, targetNucleobases, reactionNo, pressure, step,
                   rhoI, rhoP, phi, tOF, origStdOut, waterConc=None):
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

    # Calculate amounts of nucleobase(s)
    amounts = []
    iter_temps_amounts.counter = 0
    with tqdm.tqdm(total=sum([len(nos) for nos in reactionNos]),
                   desc=' 2 of 2',
                   leave=False, position=0,
                   file=origStdOut, dynamic_ncols=True) as pbar:
        for i in range(len(targetNucleobases)):
            amountsNucleobase = []
            for j in range(len(reactionNos[i])):
                amountsNucleobaseNo = iter_temps_amounts(
                                          targetNucleobases[i],
                                          reactionNos[i][j],
                                          pressure,
                                          temps,
                                          recursion=True, numbering=True,
                                          waterConc=waterConc)

#                if targetNucleobases[i] == 'ribose':
#                    amountsNucleobaseNo *= yieldRibose

                # Unit conversion
                # from [mol nucleobase/mol water]
                # to [kg nucleobase/kg planetesimal] in [ppb]
                amountsNucleobaseNo *= ((GibbsEnergy.molar_mass(
                                             targetNucleobases[i])
                                         / GibbsEnergy.molar_mass('H2O'))
                                        * (rhoI / rhoP)
                                        * phi
                                        * 1e9)

                if not targetNucleobases[i] == 'ribose':
                    amountsNucleobase.append(amountsNucleobaseNo)
                else:
                    for catalyst in range(len(yieldRibose)):
                        amountsNucleobase.append(amountsNucleobaseNo
                                                 * yieldRibose[catalyst][1])

                pbar.update(1)
            amounts.append(amountsNucleobase)

    # Close log files
    ChemApp.close_file(err)

    # Plotting
    # Styling
    plt.rcParams["figure.figsize"] = (7.101*0.5, 7.101*0.5/1.61803)
    linestyles = [(0, (2, 10)), (2, (2, 10)), (4, (2, 10)), (6, (2, 10)),
                  (8, (2, 10)), (10, (2, 10))]
    colors = plt.cm.viridis(np.linspace(1/6, 5/6, 5))
    colors = np.concatenate(([[0.01, 0.01, 0.01, 1.]], colors))
    markers = ['x', 'd', 'o', '^', 's', 'p']
    markeverys = [(0, 90), (15, 90), (30, 90), (45, 90), (60, 90), (75, 90)]
    legendHandles=[]
    styleIndex = 0

    # Plot amounts
    fig, ax1 = plt.subplots()
    for i in range(len(targetNucleobases)):
        if not targetNucleobases[i] == 'ribose':
            for j in range(len(reactionNos[i])):
                ax1.plot(time, amounts[i][j],
                         color=colors[styleIndex % len(colors)],
                         linestyle=linestyles[styleIndex % len(linestyles)],
                         marker=markers[styleIndex % len(markers)],
                         markevery=markeverys[styleIndex % len(markeverys)],
                         markersize=5)
                legendHandles.append(mlines.Line2D([], [],
                                     color=colors[styleIndex % len(colors)],
                                     marker=markers[styleIndex % len(markers)],
                                     markersize=5,
                                     label=plotLabelIndex[reactionNos[i][j]]))
                styleIndex += 1
        else:
            for j in range(len(yieldRibose)):
                ax1.plot(time, amounts[i][j],
                         color=colors[styleIndex % len(colors)],
                         linestyle=linestyles[styleIndex % len(linestyles)],
                         marker=markers[styleIndex % len(markers)],
                         markevery=markeverys[styleIndex % len(markeverys)])
                legendHandles.append(mlines.Line2D([], [],
                                     color=colors[styleIndex % len(colors)],
                                     marker=markers[styleIndex % len(markers)],
                                     label=plotLabelIndex[yieldRibose[j][0]]))
                styleIndex += 1

    # Shade regions representing experimentally found abundances in
    # C chondrites for comparision
    comp = compValues(targetNucleobases)
    for elem in comp:
        ax1.axhspan(elem[0], elem[1], alpha=0.3, color=plt.cm.viridis(0.99))
        legendHandles.append(mpatches.Patch(label=elem[2],
                                            color=plt.cm.viridis(0.99),
                                            alpha=0.3))

    ax1.set_xscale('log')
    ax1.set_xlabel('Time after formation [Myr]')#, fontsize=fontSize)
    ax1.set_xlim(np.min(time), np.max(time))#1e-3, 5e3)
    ax1.set_yscale('log')
    _, ylimTop = ax1.get_ylim()
    ax1.set_ylim(1e0, 6e2)
    ax1.set_ylabel('Molecular abundance [ppb]')#, fontsize=fontSize)
    ax1.tick_params(axis='both', which='both', top=True, direction='in')#,
#                    labelsize=fontSize)

    # Plot temperatures
    color = plt.cm.viridis(0)
    ax2 = ax1.twinx()
    ax2.plot(time, temps, color=color)
    legendHandles.append(mlines.Line2D([], [],
                         color=color,
                         label=r'Temperature $T_{\mathrm{core}}$'))
    #_, ylimTop = ax2.get_ylim()
    ax2.set_ylim(200, 525)
    ax2.set_ylabel(r'$T_{\mathrm{core}}$ [K]', color=color)#, fontsize=fontSize)
    ax2.tick_params(axis='y', direction='in', labelcolor=color)#,
#                    labelsize=fontSize)
    ax2.set_zorder(ax1.get_zorder() + 1)
    ax1.patch.set_visible(False)

    radii = tempsData[1:, 0] * 1e-3
#    titleStr = ('Temporal evolution of molecular abundances in the center of '
#                + 'the planetesimal.\nProperties of planetesimal: '
#    titleStr = ('Properties of planetesimal: '
#                + r'Radius = '
#                + '{}'.format(int(math.ceil(np.max(radii))))
#                + r' km, $\rho_{\mathrm{rock}}$ = '
#                + str(rhoP)
#                + r' kg/m$^3$, $\rho_{\mathrm{ice}}$ = '
#                + str(rhoI)
#                + r' kg/m$^3$, porosity = '
#                + str(phi)
#                + ',\ntime of formation after CAI = '
#                + str(tOF)
#                + ' Myr.')
#    if 'ribose' in targetNucleobases:
#        titleStr += ('\nExperimentally found yield of ribose within 5C sugars '
#                     + 'used: '
#                     + str(yieldRibose)
#                     + '.')
#    if not waterConc == None:
#        titleStr += ('\nInitial water concentration changed by factor of '
#                     + str(waterConc)
#                     + '.')
#    plt.title(titleStr, fontdict = {'fontsize': fontSize})

    nucleobaseTextStr = ''
    if len(targetNucleobases) == 1:
        nucleobaseTextStr = (token(targetNucleobases[0])
                             + ' Synthesis')
    elif len(targetNucleobases) == 2:
        nucleobaseTextStr = (token(targetNucleobases[0])
                             + ' and '
                             + token(targetNucleobases[1])
                             + ' Synthesis')
    else:
        for nucleobase in targetNucleobases[:-1]:
            nucleobaseTextStr += (token(nucleobase) + ', ')
        nucleobaseTextStr = nucleobaseTextStr[:-2]
        nucleobaseTextStr += (' and '
                              + token(targetNucleobases[-1])
                              + ' Synthesis')

    ax1.text(0.025, 0.915, r'\textbf{' + nucleobaseTextStr + '}',
             transform=ax1.transAxes)#, fontsize=fontSize)#, fontweight='bold'

    # Legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        ax2.legend(lines1 + lines2, labels1 + labels2, handles=legendHandles,
                   loc=3, prop={'size': 5},#-3},
                   bbox_to_anchor=(0.01, 0.01,  0.98, 0.86))

    fig.tight_layout(pad=0)

    # Save plot
    plotDir = pathlib.Path('.') / 'results'
    plotDir.mkdir(exist_ok=True)
    nucleobaseFileNameStr = ''
    for nucleobase in targetNucleobases:
        nucleobaseFileNameStr += nucleobase + '_'
    plotPath = plotDir / (nucleobaseFileNameStr
                          + str(int(pressure))
                          + 'bar_time_iter_amounts_radius_'
                          + '{}'.format(int(math.ceil(np.max(radii))))
                          + 'km.pdf')

    plt.savefig(str(plotPath),# bbox_inches='tight', pad_inches=0,
                transparent=True)


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


def plot_peak_temps_time_iter(tempsData, targetNucleobases, reactionNos,
                              pressure, step, rhoI, rhoP, phi, tOF):

    ####
    # Radial distribution
    ####

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
    for i in range(len(targetNucleobases)):
        amountsNucleobase = []
        for j in range(len(reactionNos[i])):
            amountsNucleobaseNo = iter_temps_amounts(
                                      targetNucleobases[i],
                                      reactionNos[i][j],
                                      pressure,
                                      temps,
                                      waterConc=waterConc)

            # Unit conversion
            # from [mol nucleobase/mol water]
            # to [kg nucleobase/kg planetesimal] in [ppb]
            amountsNucleobaseNo *= ((GibbsEnergy.molar_mass(
                                         targetNucleobases[i])
                                     / GibbsEnergy.molar_mass('H2O'))
                                    * (rhoI / rhoP)
                                    * phi
                                    * 1e9)

            if not targetNucleobases[i] == 'ribose':
                amountsNucleobase.append(amountsNucleobaseNo)
            else:
                for catalyst in range(len(yieldRibose)):
                    amountsNucleobase.append(amountsNucleobaseNo
                                             * yieldRibose[catalyst][1])

        amounts.append(amountsNucleobase)

    # Close log files
    ChemApp.close_file(err)

    # Plotting
    # Styling
    plt.rcParams["figure.figsize"] = (7.101, 2.7)
    #plt.rcParams["figure.figsize"] = (7.26913, 2.8)
    linestyles = [(0, (2, 10)), (2, (2, 10)), (4, (2, 10)), (6, (2, 10)),
                  (8, (2, 10)), (10, (2, 10))]
    #linestyles = [(0, (2, 8)), (2, (2, 8)), (4, (2, 8)), (6, (2, 8))]
    colors = plt.cm.viridis(np.linspace(0, 5/6, 6))
    #colors = plt.cm.viridis(np.linspace(0, 3/4, 4))
    color_temp = colors[2]
    #color_temp = colors[1]
    colors = np.delete(colors, 2, 0)
    #colors = np.delete(colors, 1, 0)
    colors = np.concatenate(([[0.01, 0.01, 0.01, 1.]], colors))
    markers = ['x', 'd', 'o', '^', 's', 'p']
    markeverys = [(0, 18), (3, 18), (6, 18), (9, 18), (12, 18), (15, 18)]
    #markeverys = [(0, 12), (3, 12), (6, 12), (9, 12)]#, (12, 18), (15, 18)]
    legendHandles=[]
    styleIndex = 0


    # Plot amounts
    fig, (ax1, ax3) = plt.subplots(ncols=2, sharey=True)
    for i in range(len(targetNucleobases)):
        if not targetNucleobases[i] == 'ribose':
            for j in range(len(reactionNos[i])):
                ax1.plot(radii, amounts[i][j],
                         color=colors[styleIndex % len(colors)],
                         linestyle=linestyles[styleIndex % len(linestyles)],
                         marker=markers[styleIndex % len(markers)],
                         markevery=markeverys[styleIndex % len(markeverys)],
                         markersize=5)
                legendHandles.append(mlines.Line2D([], [],
                                     color=colors[styleIndex % len(colors)],
                                     linestyle=(0, (2.1, 5.1)),
                                     marker=markers[styleIndex % len(markers)],
                                     markersize=5,
                                     label=plotLabelIndex[reactionNos[i][j]]))
                styleIndex += 1
        else:
            for j in range(len(yieldRibose)):
                ax1.plot(radii, amounts[i][j],
                         color=colors[styleIndex % len(colors)],
                         linestyle=linestyles[styleIndex % len(linestyles)],
                         marker=markers[styleIndex % len(markers)],
                         markevery=markeverys[styleIndex % len(markeverys)],
                         markersize=5)
                legendHandles.append(mlines.Line2D([], [],
                                     color=colors[styleIndex % len(colors)],
                                     linestyle=(0, (2.1, 5.1)),
                                     marker=markers[styleIndex % len(markers)],
                                     markersize=5,
                                     label=plotLabelIndex[yieldRibose[j][0]]))
                styleIndex += 1

    # Shade regions representing experimentally found abundances in
    # C chondrites for comparison
    comp = compValues(targetNucleobases)
    for elem in comp:
        ax1.axhspan(elem[0], elem[1], alpha=0.3, color=plt.cm.viridis(0.99))
        legendHandles.append(mpatches.Patch(label=elem[2],
                                            color=plt.cm.viridis(0.99),
                                            alpha=0.3))

    print(math.ceil(np.max(radii)))
    ax1.set_xlim(0, math.ceil(np.max(radii)))
    ax1.set_xlabel('Radius [km]')
    ax1.set_yscale('log')
    #ax1.set_ylim(1e0, 6e2)
    ax1.set_ylim(1e0, 6e2)
    ax1.set_ylabel('Molecular abundance [ppb]')
    ax1.tick_params(axis='both', which='both', top=True, direction='in')


    # Plot temperatures
    color = color_temp
    ax2 = ax1.twinx()
    ax2.plot(radii, temps, color=color)
    legendHandles.append(mlines.Line2D([], [],
                         color=color,
                         label=r'$T_{\mathrm{max}}$'))
    #ax2.set_xlim(0,np.max(radii))
    ax2.set_ylim(200, 525)
    ax2.tick_params(axis='y', direction='in', color=color)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax2.set_zorder(ax1.get_zorder() + 1)
    ax1.patch.set_visible(False)


    # Create text label to give overview over targeted molecules
    nucleobaseTextStr = ''
    if len(targetNucleobases) == 1:
        nucleobaseTextStr = (token(targetNucleobases[0])
                             + ' Synthesis')
    elif len(targetNucleobases) == 2:
        nucleobaseTextStr = (token(targetNucleobases[0])
                             + ' and '
                             + token(targetNucleobases[1])
                             + ' Synthesis')
    else:
        for nucleobase in targetNucleobases[:-1]:
            nucleobaseTextStr += (token(nucleobase) + ', ')
        nucleobaseTextStr = nucleobaseTextStr[:-2]
        nucleobaseTextStr += (' and '
                              + token(targetNucleobases[-1])
                              + ' Synthesis')

    ax1.text(0.025, 0.915, r'\textbf{' + nucleobaseTextStr + '}',
             transform=ax1.transAxes)



    ####
    # Time iteration
    ####

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

    # Calculate amounts of nucleobase(s)
    amounts = []
    iter_temps_amounts.counter = 0
    for i in range(len(targetNucleobases)):
        amountsNucleobase = []
        for j in range(len(reactionNos[i])):
            amountsNucleobaseNo = iter_temps_amounts(
                                      targetNucleobases[i],
                                      reactionNos[i][j],
                                      pressure,
                                      temps,
                                      recursion=True, numbering=True,
                                      waterConc=waterConc)

            # Unit conversion
            # from [mol nucleobase/mol water]
            # to [kg nucleobase/kg planetesimal] in [ppb]
            amountsNucleobaseNo *= ((GibbsEnergy.molar_mass(
                                         targetNucleobases[i])
                                     / GibbsEnergy.molar_mass('H2O'))
                                    * (rhoI / rhoP)
                                    * phi
                                    * 1e9)

            if not targetNucleobases[i] == 'ribose':
                amountsNucleobase.append(amountsNucleobaseNo)
            else:
                for catalyst in range(len(yieldRibose)):
                    amountsNucleobase.append(amountsNucleobaseNo
                                             * yieldRibose[catalyst][1])

        amounts.append(amountsNucleobase)

    # Close log files
    ChemApp.close_file(err)

    # Plotting
    # Styling
    linestyles = [(0, (2, 10)), (2, (2, 10)), (4, (2, 10)), (6, (2, 10)),
                  (8, (2, 10)), (10, (2, 10))]
    #linestyles = [(0, (2, 8)), (2, (2, 8)), (4, (2, 8)), (6, (2, 8))]
    colors = plt.cm.viridis(np.linspace(0, 5/6, 6))
    #colors = plt.cm.viridis(np.linspace(0, 3/4, 4))
    color_temp = colors[2]
    #color_temp = colors[1]
    colors = np.delete(colors, 2, 0)
    #colors = np.delete(colors, 1, 0)
    colors = np.concatenate(([[0.01, 0.01, 0.01, 1.]], colors))
    markers = ['x', 'd', 'o', '^', 's', 'p']
    markeverys = [(0, 90), (15, 90), (30, 90), (45, 90), (60, 90), (75, 90)]
    #markeverys = [(0, 56), (14, 56), (28, 56), (42, 56)]#, (60, 90), (75, 90)]
    styleIndex = 0

    # Plot amounts
    for i in range(len(targetNucleobases)):
        if not targetNucleobases[i] == 'ribose':
            for j in range(len(reactionNos[i])):
                ax3.plot(time, amounts[i][j],
                         color=colors[styleIndex % len(colors)],
                         linestyle=linestyles[styleIndex % len(linestyles)],
                         marker=markers[styleIndex % len(markers)],
                         markevery=markeverys[styleIndex % len(markeverys)],
                         markersize=5)
                #legendHandles.append(mlines.Line2D([], [],
                #                     color=colors[styleIndex % len(colors)],
                #                     marker=markers[styleIndex % len(markers)],
                #                     markersize=5,
                #                     label=plotLabelIndex[reactionNos[i][j]]))
                styleIndex += 1
        else:
            for j in range(len(yieldRibose)):
                ax3.plot(time, amounts[i][j],
                         color=colors[styleIndex % len(colors)],
                         linestyle=linestyles[styleIndex % len(linestyles)],
                         marker=markers[styleIndex % len(markers)],
                         markevery=markeverys[styleIndex % len(markeverys)],
                         markersize=5)
                #legendHandles.append(mlines.Line2D([], [],
                #                     color=colors[styleIndex % len(colors)],
                #                     marker=markers[styleIndex % len(markers)],
                #                     label=plotLabelIndex[yieldRibose[j][0]]))
                styleIndex += 1

    # Shade regions representing experimentally found abundances in
    # C chondrites for comparison
    comp = compValues(targetNucleobases)
    for elem in comp:
        ax3.axhspan(elem[0], elem[1], alpha=0.3, color=plt.cm.viridis(0.99))
        #legendHandles.append(mpatches.Patch(label=elem[2],
        #                                    color=plt.cm.viridis(0.99),
        #                                    alpha=0.3))

    ax3.set_xscale('log')
    ax3.set_xlabel('Time after formation [Myr]')
    ax3.set_xlim(np.min(time), np.max(time))

    ax3.tick_params(axis='x', which='minor', labelbottom=False)

    ax3.set_yscale('log')
    ax3.set_ylim(1e0, 6e2)
    #ax1.set_ylim(1e0, 8e1)
    #ax1.set_ylim(1e0, 1e4)
    #ax3.set_ylabel('Molecular abundance [ppb]')
    ax3.tick_params(axis='both', which='both', top=True, direction='in')

    # Plot temperatures
    color = color_temp
    ax4 = ax3.twinx()
    ax4.plot(time, temps, color=color, ls=(0, (1, 1)))
    legendHandles.append(mlines.Line2D([], [],
                         color=color, ls=(0, (1, 1)),
                         label=r'$T_{\mathrm{core}}$'))
    ax4.set_ylim(200, 525)
    ax4.set_ylabel(r'$T$ [K]', color=color)
    ax4.tick_params(axis='y', direction='in', color=color, labelcolor=color)
    ax4.set_zorder(ax3.get_zorder() + 1)
    ax3.patch.set_visible(False)

    radii = tempsData[1:, 0] * 1e-3

    #nucleobaseTextStr = ''
    #if len(targetNucleobases) == 1:
    #    nucleobaseTextStr = (token(targetNucleobases[0])
    #                         + ' Synthesis')
    #elif len(targetNucleobases) == 2:
    #    nucleobaseTextStr = (token(targetNucleobases[0])
    #                         + ' and '
    #                         + token(targetNucleobases[1])
    #                         + ' Synthesis')
    #else:
    #    for nucleobase in targetNucleobases[:-1]:
    #        nucleobaseTextStr += (token(nucleobase) + ', ')
    #    nucleobaseTextStr = nucleobaseTextStr[:-2]
    #    nucleobaseTextStr += (' and '
    #                          + token(targetNucleobases[-1])
    #                          + ' Synthesis')

    ax1.text(0.915, 0.915, '(a)', fontweight='bold',
             transform=ax1.transAxes)
    ax3.text(0.025, 0.915, '(b)',
             transform=ax3.transAxes)

    # Legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax3.get_legend_handles_labels()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        ax2.legend(lines1 + lines2, labels1 + labels2, handles=legendHandles,
                   #loc=(0.19, 0.027),
                   #loc=(0.022, 0.027),
                   loc=(0.4, 0.027),
                   #loc='best',
                   prop={'size': 7},
                   borderpad=0.3,
                   handletextpad=0.4)#,
                   #bbox_to_anchor=(0.022, 0.027,  0.973, 0.86))

    fig.tight_layout(pad=0)
    plt.subplots_adjust(wspace=0)
    #fig.subplots_adjust(top=0.975)

    # Save plot
    plotDir = pathlib.Path('.') / 'results'
    plotDir.mkdir(exist_ok=True)
    nucleobaseFileNameStr = ''
    for nucleobase in targetNucleobases:
        nucleobaseFileNameStr += nucleobase + '_'
    plotPath = plotDir / (nucleobaseFileNameStr
                          + str(int(pressure))
                          + 'bar_peak_temps_time_iter_amounts_radius_'
                          + '{}'.format(int(math.ceil(np.max(radii))))
                          + 'km.pdf')

    plt.savefig(str(plotPath), transparent=True)



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
# Dictionaries providing reaction numbers and legend labels for plots
####

# Dictionary with all possible reaction numbers
reactionIndex = {'adenine'  : [1, 3, 4, 6, 7, 8],
                 'uracil'   : [29, 32],
                 'cytosine' : [43, 44],
                 'guanine'  : [51, 54],
                 'thymine'  : [58, 62],
                 'ribose'   : [101],#100,
                 'C128'     : [200],
                 'adenosine': [110],
                 'uridine'  : [111],
                 'cytidine' : [112],
                 'guanosine': [113],
                 'thymidine': [114]}

# Dictionary with plot labels for all reaction numbers
plotLabelIndex = {1: r'1.\,A\,(FT)',
                    # r'1. CO + H$_2$ + NH$_3$ $\longrightarrow$ Adenine '
                    # + '+ H$_2$O',
                  3: r'3.\,A\,(NC)',
                    #r'3. 5HCN$_{(\mathrm{aq})}$ $\longrightarrow$ Adenine$_{(\mathrm{aq})}$',
                  4: r'4.\,A\,(NC)',
                    #r'4. HCN + NH$_3$ $\longrightarrow$ Adenine',
                  6: r'6.\,A\,(NC)',
                    #r'6. 5CO + 5NH$_3$ $\longrightarrow$ Adenine '
                    # + '+ 5H$_2$O',
                  7: r'7.\,A\,(NC)',
                    #r'7. HCN + H$_2$O $\longrightarrow$ Adenine',
                  8: r'8.\,A\,(NC)',
                    #r'8. HCN + NH$_3$ + H$_2$O $\longrightarrow$ Adenine',
                  29: r'29.\,U\hspace{17.7pt}(NC)',
                     #r'29. 2HCN$_{(\mathrm{aq})}$ + 2CH$_2$O$_{(\mathrm{aq})}$ '
                     # + r'$\longrightarrow$ Uracil$_{(\mathrm{aq})}$ + H$_{2(\mathrm{aq})}$',
                  32: r'32.\,C\,$\rightarrow$\,U\,(NC)',
                      #r'32. Cytosine + H$_2$O $\longrightarrow$ Uracil + '
                      #+ r'NH$_3$',
                  43: r'43.\,C\hspace{17.9pt}(FT)',
                      #r'43. CO + H$_2$ + NH$_3$ $\longrightarrow$ Cytosine + '
                      #+ r'H$_2$O',
                  44: r'44.\,C\hspace{17.9pt}(NC)',
                      #r'44. 3HCN$_{(\mathrm{aq})}$ + CH$_2$O$_{(\mathrm{aq})}$ $\longrightarrow$'
                      #+ r' Cytosine$_{(\mathrm{aq})}$',
                  51: r'51.\,G\hspace{17.4pt}(FT)',
                      #r'51. CO + H$_2$ + NH$_3$ $\longrightarrow$ Guanine + '
                      #+ r'H$_2$O',
                  54: r'54.\,G\hspace{17.4pt}(NC)',
                      #r'54. 5HCN$_{(\mathrm{aq})}$ + H$_2$O $\longrightarrow$ '
                      #+ r'Guanine$_{(\mathrm{aq})}$ + H$_2(\mathrm{aq})$',
                  58: r'58.\,T\hspace{17.9pt}(NC)',
                      #r'58. 2HCN$_{(\mathrm{aq})}$ + 3CH$_2$O$_{(\mathrm{aq})}$ '
                      #+ r'$\longrightarrow$ Thymine$_{(\mathrm{aq})}$ + H$_2$O',
                  62: r'62.\,U\,$\rightarrow$\,T\,(NC)',
                      #r'62. Uracil + CH$_2$O + CH$_2$O$_2$ + H$_2$O '
                      #+ r'$\longrightarrow$ Thymine',
                  100: r'100. 5CH$_2$O$_{(\mathrm{aq})}$ $\longrightarrow$ '
                       + r'Ribose$_{(\mathrm{aq})}$',
                  101: r'101. $\ce{CH2O_{(\mathrm{aq})} + 2C2H4O2_{(\mathrm{aq})} -> \text{Ribose}_{(\mathrm{aq})}}$',
                  'Ca(OH)2': r'101a. $\ce{CH2O_{(\mathrm{aq})} + 2C2H4O2_{(\mathrm{aq})} ->[Ca(OH)2] C5H10O5_{(\mathrm{aq})}}$',
                  'CaCO3': r'101b. $\ce{CH2O_{(\mathrm{aq})} + 2C2H4O2_{(\mathrm{aq})} ->[CaCO3\hspace{4.5pt}] C5H10O5_{(\mathrm{aq})}}$',
                  'K2CO3': r'101c. $\ce{CH2O_{(\mathrm{aq})} + 2C2H4O2_{(\mathrm{aq})} ->[K2CO3\hspace{4.7pt}] C5H10O5_{(\mathrm{aq})}}$',
                  'KOH': r'101d. $\ce{CH2O_{(\mathrm{aq})} + 2C2H4O2_{(\mathrm{aq})} ->[KOH\hspace{8.9pt}] C5H10O5_{(\mathrm{aq})}}$',
                  110: r'110. Adenine$_{(\mathrm{aq})}$ + Ribose$_{(\mathrm{aq})}$ '
                       + r'$\longrightarrow$ Adenosine$_{(\mathrm{aq})}$ + H$_2$O',
                  111: r'111. Uracil$_{(\mathrm{aq})}$ + Ribose$_{(\mathrm{aq})}$ '
                       + r'$\longrightarrow$ Uridine$_{(\mathrm{aq})}$ + H$_2$O',
                  112: r'112. Cytosine$_{(\mathrm{aq})}$ + Ribose$_{(\mathrm{aq})}$ '
                       + r'$\longrightarrow$ Cytidine$_{(\mathrm{aq})}$ + H$_2$O',
                  113: r'113. Guanine$_{(\mathrm{aq})}$ + Ribose$_{(\mathrm{aq})}$ '
                       + r'$\longrightarrow$ Guanosine$_{(\mathrm{aq})}$ + H$_2$O',
                  114: r'114. Thymine$_{(\mathrm{aq})}$ + Ribose$_{(\mathrm{aq})}$ '
                       + r'$\longrightarrow$ Thymidine$_{(\mathrm{aq})}$ + H$_2$O',
                  200: r'200. 128CH$_2$O$_{(\mathrm{aq})}$ $\longrightarrow$ '
                       + r'C$_{128}$H$_{68}$O$_{7(\mathrm{s})}$ + 54OH$^-_{(\mathrm{aq})}$ + '
                       + r'67H$_2$O'}



####
# --------------------------- Main Program ------------------------------------
####

####
# Command line arguments
####

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
# Define properties of asteroid
####

rhoI = 917  # density of water ice [kg/m^3]
rhoP = 3000 # density of planetesimal [kg/m^3]
phi  = 0.2  # porosity
#tOF  = 2.5  # time of formation after CAI [Myr]

# Use to artificially change initial water concentration
# If not desired, set to None
waterConc = None#0.0005

# Set experimentally found yield of ribose within 5C sugars
#yieldRibose = 0.12#0.07
yieldRibose = [['Ca(OH)2', 0.04096157],
               ['CaCO3'  , 0.03446486],
               ['K2CO3'  , 0.027312665],
               ['KOH'    , 0.024244693]]

####
# Read temperature input file(s)
####

tempsDir  = pathlib.Path('.') / 'temps_input'
if not tempsDir.is_dir():
    errorStr =  'Directory \'' + str(tempsDir) + '\'necessary '
    errorStr += 'to read input temperature file(s) from not found'
    raise NotADirectoryError(errorStr)
tempsPaths = tempsDir.glob('**/*')
tempsFiles = [x for x in tempsPaths]

with std_out_err_redirect_tqdm() as origStdOut:
    for i in tqdm.trange(len(tempsFiles), desc='overall', leave=False,
                         position=1,
                         file=origStdOut, dynamic_ncols=True):
        tempsData = np.genfromtxt(tempsFiles[i], delimiter=',', unpack=True)

        # Find time of formation after CAI [Myr] from file name
        tOF = tempsFiles[i].name
        tOF = float(tOF[tOF.find('km_') + 3:tOF.find('Myr.csv')].replace('_',
                                                                       '.'))
        ####
        # Calculate amounts and plot
        ####

        # Disable warning of matplotlib for too many open figures
        plt.rcParams.update({'figure.max_open_warning': 0})

        #plot_peak_temps(tempsData, targetNucleobases, reactionNos, pressure,
        #                rhoI, rhoP, phi, tOF, origStdOut, waterConc)

        step = 1
        #plot_time_iter(tempsData, targetNucleobases, reactionNos, pressure,
        #               step,
        #               rhoI, rhoP, phi, tOF, origStdOut, waterConc)
        plot_peak_temps_time_iter(tempsData, targetNucleobases, reactionNos,
                                  pressure, step, rhoI, rhoP, phi, tOF)

print('\n*** DONE ***')
