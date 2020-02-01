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

#include "cacint.h"


// If a ChemApp error occures, this function reports the error number and the
// routine it occured in before it exits the program.
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
    // variable declarations for use with ChemApp routines
    ////

    // integer variables
    // error return variable
    LI noErr,
    // version number of ChemApp
    caVers,
    // used to return information from 'tqsiz' about capabilities of used
    // ChemApp version
    la, lb, lc, lf, lh, li, lj, ll, lm, ln, lo, lp, lq, lr,
    // FORTRAN unit number under which the thermochemical data-files will be
    // opened
    unitNo,
    // user ID, ChemApp license holder and dongle info
    hasPid, edMon, edYear;

    // floating point variables


    // string variables
    // user ID, ChemApp license holder and dongle info
    char id[255], name[80], pid[TQSTRLEN], hasPt[TQSTRLEN];


    ////
    // ChemApp starter routines
    ////

    // initialize ChemApp
    tqini(&noErr);
    if (noErr) abortProg(__LINE__, "tqini", noErr);

    // print ChemApp copyright message
    tqcprt(&noErr);
    if (noErr) abortProg(__LINE__, "tqcprt", noErr);
    tqerr((CHP)TQERRMSG, &noErr);
    if (noErr) abortProg(__LINE__, "tqerr", noErr);
    for (std::size_t i(0); i < 3; ++i)
        std::cout << TQERRMSG[i] << std::endl;
    std::cout << std::endl;

    // get ChemApp version number
    tqvers(&caVers, &noErr);
    if (noErr) abortProg(__LINE__, "tqvers", noErr);
    std::cout << "ChemApp version is: " << caVers << "\n\n";

    // get the licensee's user ID, name and program ID
    tqgtid(id, &noErr);
    if (noErr) abortProg(__LINE__, "tqgtid", noErr);
    tqgtnm(name, &noErr);
    if (noErr) abortProg(__LINE__, "tqgtnm", noErr);
    tqgtpi(pid, &noErr);
    if (noErr) abortProg(__LINE__, "tqgtpi", noErr);

    // print the licensee's user ID, name and program ID
    std::cout << "Licensee\'s user ID: " << id
              << "\nLicensee\'s name: " << name
              << "\nProgram ID: " << pid << "\n\n";

    // The following pieces of information are only meaningful if a version
    // of ChemApp is used that requires a dongle (hardware key).
    // get the HASP dongle type and id
    tqgthi(hasPt, &hasPid, &noErr);
    if (noErr) abortProg(__LINE__, "tqgthi", noErr);
    // get the ChemApp license expiration date (month and year)
    tqgted(&edMon, &edYear, &noErr);
    if (noErr) abortProg(__LINE__, "tqgted", noErr);
    // print info if HASP dongle is used
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


    // get array sizes to get information about capabilities of ChemApp version
    // used
    tqnsiz(&la, &lb, &lc, &lf, &lh, &li, &lj, &ll, &lm, &ln, &lo, &lp, &lq,
           &lr, &noErr);
    if (noErr) abortProg(__LINE__, "tqnsiz", noErr);
    std::cout << "Maximum number of constituents:                    "   << la
              << "\nMaximum number of system components:               " << lb
              << "\nMaximum number of mixture phases:                  " << lc
              << "\nMaximum number of sublattices for a mixture phase: " << lf
              << "\n\n";

    // get value for input/output option
    // determine which FORTRAN unit is used by default by 'tqrfil' for reading
    // the thermochemical data-file, and use the returned number to
    // subsequently open and close the data-file
    tqgio((char*)"FILE", &unitNo, &noErr);
    if (noErr) abortProg(__LINE__, "tqgio", noErr);
    std::cout << "The thermochemical data will be read from the file "
              << "associated with unit " << unitNo << "\n\n";

    ////
    // ChemApp data collection
    ////



    return 0;
}
