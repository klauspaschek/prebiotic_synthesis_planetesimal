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
 */

#include <iostream>
#include <string>
#include <cstddef>
#include <cstdlib>
#include <iomanip>

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
    char fileName[32] = "moleculeData100barAdenine-3.dat";
    tqopna(fileName, unitNo, &err);
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
    for (int i(1); i <= nSCom; ++i)
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
    for (int i(1); i <= nPhase; ++i)
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
    for (int i(1); i <= nPhase; ++i)
    {
        tqgnp(i, pName, &err);
        if (err) abortProg(__LINE__, "tqgnp", err);
        tqnopc(i, &nPCon, &err);
        if (err) abortProg(__LINE__, "tqnopc", err);
        for (int j(1); j <= nPCon; ++j)
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
    for (int i(1); i <= nPhase; ++i)
    {
        tqgnp(i, pName, &err);
        if (err) abortProg(__LINE__, "tqgnp", err);
        tqnopc(i, &nPCon, &err);
        if (err) abortProg(__LINE__, "tqnopc", err);
        // We only need to check mixture phases
        if (nPCon > 1)
        {
            for (int j(1); j <= nPCon; ++j)
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

    // If there were phase constituents which cannot be used as incoming
    // species, print a note
    if (nPerm > 0)
        std::cout << "Note: " << nPerm << " phase constituent(s) is/are not "
                  << "permitted as incoming species" << std::endl;


    return 0;
}
