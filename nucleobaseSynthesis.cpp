/* Strecker Synthesis ChemApp Simulation Program
 * This program makes calls to ChemApp subroutines to simulate the synthesis of
 * nucleobases for an input temperature & pressure and abundances of initial
 * molecules.
 *
 * Author:      Ben Pearce
 * Port to C++: Klaus Paschek
 * Sources:     ChemApp Programmer's Manual
 * (https://gtt-technologies.de/software/chemapp/documentation/online-manual/,
 * January 2020)
 *
 * Note: Be aware that often loop indices follow the FORTRAN index convention
 * by starting with 1, not 0, to make the ChemApp routines work properly. When
 * using them for element access within C/C++ containers, they often have to be
 * reduced by 1 so that again these containters behave correctly.
 * Example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for (LI i(1); i <= size; ++i)    // LI is int usable with ChemApp routines
 *     tq...(..., i, ..., &err);    // ChemApp routine
 *     someVariable = someContainer[i - 1];
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include <iostream>
#include <string>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iterator>
#include <utility>
#include <cstring>

// ChemApp
#include "cacint.h"


// If a ChemApp error occures, this function reports the error number and the
// routine it occured in before it exits the program
void abortProg(int lineNo, std::string funcName, LI errorNo)
{
    std::cout << "ChemApp error no. " << errorNo << " occured when calling "
              << funcName << ".\nAborting on line " << lineNo << " of \'"
              << __FILE__ << "\'." << std::endl;
    exit(errorNo);
}


int main()
{
    ////
    // Wanted nucleobase
    ////
    char nucleobaseStr[TQSTRLEN] = "C5H5N5";


    ////
    // Input/output file names
    ////

    // Input
    std::string thermoChemFileName("moleculeData100barAdenine-3.dat");
    std::string tempsFileName     ("inputPressureTemps.dat");
    std::string initConcFileName  ("initialConcentrations.dat");

    // Output
    // File containing nucleobase amounts as main result of program
    std::string amountsFileName   ("adenine100barWaterRatio-3.dat");
    // File containing molar masses and conditions for individual equilibrium
    // calculations
    std::string molMassFileName   ("moleMasses.dat");


    // Error return variable that will be checked after each call of a ChemApp
    // routine
    LI err;


    ////
    // ChemApp starter routines
    ////

    // Initialize ChemApp
    tqini(&err);
    if (err) abortProg(__LINE__, "tqini", err);

    // Print ChemApp copyright message
    tqcprt(&err);
    if (err) abortProg(__LINE__, "tqcprt", err);
    tqerr((CHP)TQERRMSG, &err);
    if (err) abortProg(__LINE__, "tqerr", err);
    for (int i(0); i < 3; ++i)
        std::cout << TQERRMSG[i] << std::endl;
    std::cout << std::endl;

    // Get ChemApp version number
    LI caVers;
    tqvers(&caVers, &err);
    if (err) abortProg(__LINE__, "tqvers", err);

    std::cout << "ChemApp version is: " << caVers << "\n\n";

    // Get the licensee's user ID, name and program ID
    char id[255], name[80], pid[TQSTRLEN];
    tqgtid(id, &err);
    if (err) abortProg(__LINE__, "tqgtid", err);
    tqgtnm(name, &err);
    if (err) abortProg(__LINE__, "tqgtnm", err);
    tqgtpi(pid, &err);
    if (err) abortProg(__LINE__, "tqgtpi", err);

    // Print the licensee's user ID, name and program ID
    std::cout << "Licensee\'s user ID: " << id
              << "\nLicensee\'s name: " << name
              << "\nProgram ID: " << pid << "\n\n";

    // The following pieces of information are only meaningful if a version
    // of ChemApp is used that requires a dongle (hardware key).
    // Get the HASP dongle type and id
    LI hasPid, edMon, edYear;
    char hasPt[TQSTRLEN];
    tqgthi(hasPt, &hasPid, &err);
    if (err) abortProg(__LINE__, "tqgthi", err);

    // Get the ChemApp license expiration date (month and year)
    tqgted(&edMon, &edYear, &err);
    if (err) abortProg(__LINE__, "tqgted", err);

    // Print info if HASP dongle is used
    if (hasPid)
    {
        std::cout << "HASP dongle type: " << hasPt
                  << "\nHASP dongle id: " << hasPid
                  << "\nChemApp license expiration date (month/year): "
                  << edMon << "/" << edYear << "\n\n";
    } else {
        std::cout << "This ChemApp version does not require a HASP hardware "
                  << "key (dongle)" << "\n\n";
    }

    // Get array sizes to get information about capabilities of ChemApp version
    // used
    LI la, lb, lc, lf, lh, li, lj, ll, lm, ln, lo, lp, lq, lr;
    tqnsiz(&la, &lb, &lc, &lf, &lh, &li, &lj, &ll, &lm, &ln, &lo, &lp, &lq,
           &lr, &err);
    if (err) abortProg(__LINE__, "tqnsiz", err);

    std::cout << std::setw(51) << std::left
              << "Maximum number of constituents:"   << la
              << std::setw(52) << std::left
              << "\nMaximum number of system components:" << lb
              << std::setw(52) << std::left
              << "\nMaximum number of mixture phases:" << lc
              << std::setw(52) << std::left
              << "\nMaximum number of sublattices for a mixture phase:" << lf
              << "\n" << std::endl;


    ////
    // ChemApp data collection
    ////

    // Get value for input/output option
    // Determine which FORTRAN unit is used by default by 'tqrfil' for reading
    // the thermochemical data-file, and use the returned number to
    // subsequently open and close the data-file
    LI unitNo;
    tqgio((char*)"FILE ", &unitNo, &err);
    if (err) abortProg(__LINE__, "tqgio", err);

    std::cout << "The thermochemical data will be read from the file "
              << "associated with unit " << unitNo << "\n\n";

    // Open ASCII thermochemical data-file for reading
    tqopna((char*)thermoChemFileName.c_str(), unitNo, &err);
    if (err) abortProg(__LINE__, "tqopna", err);

    // Read data-file
    tqrfil(&err);
    if (err) abortProg(__LINE__, "tqrfil", err);

    // Close data-file
    tqclos(unitNo, &err);
    if (err) abortProg(__LINE__, "tqclos", err);

    // Get default system units
    char dstr[TQSTRLEN];
    std::cout << "DEFAULT SYSTEM UNITS:\nQuantity    Unit" << std::endl;

    tqgsu((char*)"Pressure ", dstr, &err);
    if (err) abortProg(__LINE__, "tqgsu", err);
    std::cout << "Pressure    " << dstr << std::endl;

    tqgsu((char*)"Volume ", dstr, &err);
    if (err) abortProg(__LINE__, "tqgsu", err);
    std::cout << "Volume      " << dstr << std::endl;

    tqgsu((char*)"Temperature ", dstr, &err);
    if (err) abortProg(__LINE__, "tqgsu", err);
    std::cout << "Temperature " << dstr << std::endl;

    tqgsu((char*)"Energy ", dstr, &err);
    if (err) abortProg(__LINE__, "tqgsu", err);
    std::cout << "Energy      " << dstr << std::endl;

    tqgsu((char*)"Amount ", dstr, &err);
    if (err) abortProg(__LINE__, "tqgsu", err);
    std::cout << "Amount      " << dstr << "\n" <<  std::endl;

    // Change system unit (default of 'Amount' is mol)
    // Here we change the unit for the quantity 'Amount' to gram, mainly so
    // that when we call 'tqstsc' firther down, we get the molecular mass
    // expressed in the unit g/mol
    tqcsu((char*)"Amount ", (char*)"gram ", &err);
    if (err) abortProg(__LINE__, "tqcsu", err);
    std::cout << "Amount unit set to gram for masses below\n" << std::endl;

    // Get number of system components
    LI nSCom;
    tqnosc(&nSCom, &err);
    if (err) abortProg(__LINE__, "tqnosc", err);

    std::cout << "SYSTEM COMPONENTS\nNumber of system compoents: " << nSCom
              << std::endl;

    // Print the names of the system components and their masses
    std::cout << "No.  Name of component  Mass [g/mol]" << std::endl;
    for (LI i(1); i <= nSCom; ++i)
    {
        tqgnsc(i, dstr, &err);
        if (err) abortProg(__LINE__, "tqgnsc", err);

        DB wMass;
        DB stoi[nSCom];
        tqstsc(i, stoi, &wMass, &err);
        if (err) abortProg(__LINE__, "tqstsc", err);

        std::cout << i << "    " << std::setw(19) << std::left << dstr
                  << wMass << std::endl;
    }
    std::cout << std::endl;

    // Get number of phases
    LI nPhase;
    char pName[TQSTRLEN];
    tqnop(&nPhase, &err);
    if (err) abortProg(__LINE__, "tqnop", err);

    // Print the names of the phases and their model names
    std::cout << "PHASES\nNumber of phases: " << nPhase
              << "\nNo.  Name of phase  Model" << std::endl;

    for (LI i(1); i <= nPhase; ++i)
    {
        tqgnp(i, pName, &err);
        if (err) abortProg(__LINE__, "tqgnp", err);
        tqmodl(i, dstr, &err);
        if (err) abortProg(__LINE__, "tqmodl", err);

        std::cout << i << "    " << std::setw(15) << std::left << pName
                  << dstr << std::endl;
    }
    std::cout << "(PURE: Stoichiometric condensed phase, IDMX: Ideal mixing)"
              << "\n" << std::endl;

    // Get number of phase constituents of gas phase
    LI nPCon;
    tqnopc(1, &nPCon, &err);
    if (err) abortProg(__LINE__, "tqnopc", err);

    // Print the names of the phase constituents and their Gibbs free energies
    // of formation
    std::cout << "PHASE CONSTITUENTS\nNumber of phase constituents of the "
              << "first gas phase: " << nPCon
              << "\nPh.  Name of constituent  Gibbs [J/mol]\n";

    char pCName[TQSTRLEN];
    for (LI i(1); i <= nPhase; ++i)
    {
        tqgnp(i, pName, &err);
        if (err) abortProg(__LINE__, "tqgnp", err);
        tqnopc(i, &nPCon, &err);
        if (err) abortProg(__LINE__, "tqnopc", err);

        for (LI j(1); j <= nPCon; ++j)
        {
            tqgnpc(i, j, pCName, &err);
            if (err) abortProg(__LINE__, "tqgnpc", err);

            LI nValv;
            DB valv[25];
            tqgdat(i, j, (char*)"H ", 1, &nValv, valv, &err);
            if (err) abortProg(__LINE__, "tqgdat", err);

            DB sEnthal(valv[0]);
            tqgdat(i, j, (char*)"S ", 1, &nValv, valv, &err);
            if (err) abortProg(__LINE__, "tqgdat", err);
            DB sEntrop(valv[0]);

            // G = H - T*S
            DB sGibbs(sEnthal - (298.15 * sEntrop));
            std::cout << i << "    " << std::setw(19) << std::left << pCName
                      << "  " << sGibbs << std::endl;
        }
    }
    std::cout << std::endl;

    // Change system unit back to mol
    tqcsu((char*)"Amount ", (char*)"mol ", &err);
    if (err) abortProg(__LINE__, "tqcsu", err);
    std::cout << "Amount unit set back to mol" << "\n" << std::endl;

    // Make sure everything can be used as an incoming species
    std::cout << "Mixture phase  Phase constituent  Ok/Not permitted"
              << std::endl;

    LI nPerm(0);
    for (LI i(1); i <= nPhase; ++i)
    {
        tqgnp(i, pName, &err);
        if (err) abortProg(__LINE__, "tqgnp", err);
        tqnopc(i, &nPCon, &err);
        if (err) abortProg(__LINE__, "tqnopc", err);

        // We only need to check mixture phases
        if (nPCon > 1)
        {
            for (LI j(1); j <= nPCon; ++j)
            {
                tqgnpc(i, j, pCName, &err);
                if (err) abortProg(__LINE__, "tqgnpc", err);

                // Check whether the current constituent is permitted to be
                // used as incoming species
                LI inPCIS;
                tqpcis(i, j, &inPCIS, &err);
                if (err) abortProg(__LINE__, "tqpcis", err);

                // Print a table entry, and if it is not permitted, increase
                // the total number of such phase constituents by 1
                std::cout << std::setw(15) << std::left << pName
                          << std::setw(19) << std::left << pCName;
                if (inPCIS == 0)
                {
                    ++nPerm;
                    std::cout << "Not peritted" << std::endl;
                } else {
                    std::cout << "Ok" << std::endl;
                }
            }
        }
    }
    std::cout << std::endl;

    // If there were phase constituents which cannot be used as incoming
    // species, print a note
    if (nPerm > 0)
        std::cout << "Note: " << nPerm << " phase constituent(s) is/are not "
                  << "permitted as incoming species" << std::endl;


    ////
    // Main Program
    ////

    // Read in initial temperatures to iterate program over.
    // Store them in array 'temps'
    std::ifstream tempsFile(tempsFileName);
    std::string dummyStr;
    DB pressure;
    tempsFile >> dummyStr >> pressure;
    std::istream_iterator<DB> start(tempsFile), end;
    std::vector<DB> temps(start, end);

    // Open resulting files containing nucleobase amounts and molecular masses
    std::ofstream amountsFile(amountsFileName);
    std::ofstream molMassFile(molMassFileName);

    // Header for nucleobase amounts file
    amountsFile << "Temperature [C],Amount [kg nucleobase/kg water]"
                << std::endl;

    // Begin iterative loop over all temperatures
    for (DB temperature:temps)
    {
        // Set temperature and pressure
        LI numCon;
        tqsetc((char*)"T ", 0, 0, temperature, &numCon, &err);
        if (err) abortProg(__LINE__, "tqsetc", err);
        tqsetc((char*)"P ", 0, 0, pressure, &numCon, &err);
        if (err) abortProg(__LINE__, "tqsetc", err);

        // Obtain the inital constituent concentrations
        std::ifstream initConcFile(initConcFileName);
        std::istream_iterator<DB> start(initConcFile), end;
        std::vector<DB> concs(start, end);

        // Set these initial concentrations in ChemApp
        int indexConcs(0);
        for (LI i(1); i <= nPhase; ++i)
        {
            // Find number of constituents in phase
            tqnopc(i, &nPCon, &err);
            if (err) abortProg(__LINE__, "tqnopc", err);

            for (LI j(1); j <= nPCon; ++j)
            {
                tqsetc((char*)"ia ", i, j, concs[indexConcs], &numCon, &err);
                if (err) abortProg(__LINE__, "tqsetc", err);

                ++indexConcs;
            }
        }

        // Redirect the results (output of ChemApp routines like 'tqcel' and
        // 'tqshow') to an output file for viewing
        tqcio((char*)"LIST ", 3, &err);
        if (err) abortProg(__LINE__, "tqcio", err);
        tqopen((char*)"results.dat", 3, &err);
        if (err) abortProg(__LINE__, "tqopen", err);

        // Use 'tqwstr' to write header to the output file
        tqwstr((char*)"LIST ", (char*)"RESULTS ", &err);
        if (err) abortProg(__LINE__, "tqwstr", err);

        // Calculate and write equilibrium to the output file
        // The "  " option means there will be no target variable estimate used
        // to begin the simulator at, meaning 'estim" istn't used
        DB estim[1] = { 0 };
        tqcel((char*)" ", 0, 0, estim, &err);
        if (err) abortProg(__LINE__, "tqcel", err);

        // Change 'Amount' unit to gram for calculations of ppm
        tqcsu((char*)"Amount ", (char*)"gram ", &err);
        if (err) abortProg(__LINE__, "tqcsu", err);

        // Determine molecular masses of the phase constituents
        std::vector< std::vector<DB> > molMass;
        molMass.resize(nPhase);
        for (LI i(1); i <= nPhase; ++i)
        {
            // Find number of constituents in phase
            tqnopc(i, &nPCon, &err);
            if (err) abortProg(__LINE__, "tqnopc", err);

            DB wMass;
            DB stoi[nSCom];
            for (LI j(1); j <= nPCon; ++j)
            {
                tqstpc(i, j, stoi, &wMass, &err);
                if (err) abortProg(__LINE__, "tqstpc", err);

                molMass[i - 1].push_back(wMass);
            }
        }

        // Change 'Amount' unit back to mol
        tqcsu((char*)"Amount ", (char*)"mol ", &err);
        if (err) abortProg(__LINE__, "tqcsu", err);

        // Write temperature and pressure to molecular masses output file
        DB value;
        // Temperature
        tqgetr((char*)"T ", 0, 0, &value, &err);
        if (err) abortProg(__LINE__, "tqgetr", err);

        molMassFile << "Temperature: " << value << std::endl;

        // Save current temperature to write it later to other output file
        DB currentTemp(value);

        // Pressure
        tqgetr((char*)"P ", 0, 0, &value, &err);
        if (err) abortProg(__LINE__, "tqgetr", err);

        molMassFile << "Pressure:    " << value << std::endl;

        // Place equilibrium concentrations of nucleobases into an output file
        molMassFile << "Phase constituents molar masses [g/mol]:" << std::endl;

        // If water takes part in the reaction, collect corresponding indices
        std::vector< std::pair<LI, LI> > indicesWater;
        // Also collect indices corresponding to wanted nucleobase
        std::vector< std::pair<LI, LI> > indicesNucleobase;
        // Also collect molar mass of wanted nucleobase
        DB molMassNucleobase;
        // Furthermore collect all molar masses and write them to output file
        for (LI i(1); i <= nPhase; ++i)
        {
            tqnopc(i, &nPCon, &err);
            if (err) abortProg(__LINE__, "tqnopc", err);

            for (LI j(1); j <= nPCon; ++j)
            {
                tqgnpc(i, j, pCName, &err);
                if (err) abortProg(__LINE__, "tqgnpc", err);

                molMassFile << std::setw(10) << std::left << pCName
                            << molMass[i - 1][j - 1] << std::endl;

                // Check if constituent is water and if collect corresponding
                // indices
                char h2oStr[TQSTRLEN] = "H2O";
                if (!std::strcmp(pCName, h2oStr))
                    indicesWater.push_back(std::make_pair(i, j));

                // Check if constituent is the wanted nucleobase and collect
                // corresponding indices & molar mass
                if (!std::strcmp(pCName, nucleobaseStr))
                {
                    indicesNucleobase.push_back(std::make_pair(i, j));
                    molMassNucleobase = molMass[i - 1][j - 1];
                }
            }
        }
        molMassFile << "\nEquilibrium amounts (species, mol, ppb)\n"
                    << std::endl;

        // If part of reaction, get concentration of water
        DB waterConc(100.0);
        if (indicesWater.size() != 0)
        {
            for (std::pair<LI, LI> indices:indicesWater)
            {
                DB phaseWaterConc;
                tqgetr((char*)"A ", indices.first, indices.second,
                       &phaseWaterConc, &err);
                if (err) abortProg(__LINE__, "tqgetr", err);

                waterConc += phaseWaterConc;
            }
        }

        // Collect total adenine abundance in moles
        DB nucleobaseConc(0);
        for (std::pair<LI, LI> indices:indicesNucleobase)
        {
            DB phaseNucleobaseConc;
            tqgetr((char*)"A ", indices.first, indices.second,
                   &phaseNucleobaseConc, &err);
            if (err) abortProg(__LINE__, "tqgetr", err);

            nucleobaseConc += phaseNucleobaseConc;
        }

        // Get ratio of nucleobase X to water (kg X/kg water)
        DB nucleobaseRatio( (nucleobaseConc * molMassNucleobase)
                           /(waterConc      * 18.015));

        // Write temperature (in °C) and nucleobase amount to file (main result
        // of program)
        amountsFile << currentTemp - 273.15 << "," << nucleobaseRatio
                    << std::endl;
    }


    ////
    // Closing routines
    ////

    // Close file containing redirected results of ChemApp routines like
    // 'tqcel' and 'tqshow'
    tqclos(3, &err);
    if (err) abortProg(__LINE__, "tqclos", err);

    // Redirect the output back to unit number 6 (standard output unit)
    tqcio((char*)"LIST ", 6, &err);
    if (err) abortProg(__LINE__, "tqcio", err);


    return 0;
}
