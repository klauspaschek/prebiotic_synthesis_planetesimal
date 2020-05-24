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
 *     someVariable = someContainer[i - 1];    // C/C++ container
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
#include <cctype>
#include <algorithm>
#include <exception>
#include <locale>

// ChemApp
#include "cacint.h"

// pybind11
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
namespace py = pybind11;
using namespace pybind11::literals;



// If a ChemApp error occures, this function reports the error number and the
// routine it occured in before it exits the program
void abort_prog(int lineNo, std::string funcName, LI errorNo)
{
    std::cout << "ChemApp error no. "
              << errorNo
              << " occured when calling "
              << funcName
              << ".\nAborting on line "
              << lineNo
              << " of \'"
              << __FILE__
              << "\'.\n";
    exit(errorNo);
}


// Custom exception class
class ErrorStr : public std::exception
{
    private:
        std::string m_message;

    public:
        ErrorStr(std::string message)
        : m_message(message)
        {}
        ErrorStr() = delete;

        virtual const char* what() const noexcept override
        {
            return m_message.c_str();
        }
};


// Use comma ',' as delimiter in stream classes of std::
// For example to read '.csv' files
struct CSVDelimiter : std::ctype<char>
{
    CSVDelimiter()
    : std::ctype<char>(get_table())
    {}

    static mask const* get_table()
    {
        static mask rc[table_size];
        rc[','] = std::ctype_base::space;
        rc['\n'] = std::ctype_base::space;
        return &rc[0];
    }
};


// Startup of ChemApp including initialization and printing of copyright
// message & capabilities of linked ChemApp version
void start(LI& err)
{
    // Initialize ChemApp
    tqini(&err);
    if (err) abort_prog(__LINE__, "tqini", err);

    // Print ChemApp copyright message
    tqcprt(&err);
    if (err) abort_prog(__LINE__, "tqcprt", err);
    tqerr((CHP)TQERRMSG, &err);
    if (err) abort_prog(__LINE__, "tqerr", err);
    for (int i(0); i < 3; ++i)
        std::cout << TQERRMSG[i] << "\n";
    std::cout << "\n";

    // Get ChemApp version number
    LI caVers;
    tqvers(&caVers, &err);
    if (err) abort_prog(__LINE__, "tqvers", err);

    std::cout << "ChemApp version is: " << caVers << "\n\n";

    // Get the licensee's user ID, name and program ID
    char id[255], name[80], pid[TQSTRLEN];
    tqgtid(id, &err);
    if (err) abort_prog(__LINE__, "tqgtid", err);
    tqgtnm(name, &err);
    if (err) abort_prog(__LINE__, "tqgtnm", err);
    tqgtpi(pid, &err);
    if (err) abort_prog(__LINE__, "tqgtpi", err);

    // Print the licensee's user ID, name and program ID
    std::cout << "Licensee\'s user ID: "
              << id
              << "\nLicensee\'s name: "
              << name
              << "\nProgram ID: "
              << pid
              << "\n\n";

    // The following pieces of information are only meaningful if a version
    // of ChemApp is used that requires a dongle (hardware key).
    // Get the HASP dongle type and id
    LI hasPid, edMon, edYear;
    char hasPt[TQSTRLEN];
    tqgthi(hasPt, &hasPid, &err);
    if (err) abort_prog(__LINE__, "tqgthi", err);

    // Get the ChemApp license expiration date (month and year)
    tqgted(&edMon, &edYear, &err);
    if (err) abort_prog(__LINE__, "tqgted", err);

    // Print info if HASP dongle is used
    if (hasPid)
    {
        std::cout << "HASP dongle type: "
                  << hasPt
                  << "\nHASP dongle id: "
                  << hasPid
                  << "\nChemApp license expiration date (month/year): "
                  << edMon
                  << "/"
                  << edYear
                  << "\n\n";
    } else {
        std::cout << "This ChemApp version does not require a HASP hardware "
                  << "key (dongle)\n\n";
    }

    // Get array sizes to get information about capabilities of ChemApp version
    // used
    LI la, lb, lc, lf, lh, li, lj, ll, lm, ln, lo, lp, lq, lr;
    tqnsiz(&la, &lb, &lc, &lf, &lh, &li, &lj, &ll, &lm, &ln, &lo, &lp, &lq,
           &lr, &err);
    if (err) abort_prog(__LINE__, "tqnsiz", err);

    std::cout << std::setw(51)
              << std::left
              << "Maximum number of constituents:"
              << la
              << std::setw(52)
              << std::left
              << "\nMaximum number of system components:"
              << lb
              << std::setw(52)
              << std::left
              << "\nMaximum number of mixture phases:"
              << lc
              << std::setw(52)
              << std::left
              << "\nMaximum number of sublattices for a mixture phase:"
              << lf
              << "\n\n";
}


void read_data(std::string fileName, LI& err)
{
    // Get value for input/output option
    // Determine which FORTRAN unit is used by default by 'tqrfil' for reading
    // the thermochemical data-file, and use the returned number to
    // subsequently open and close the data-file
    LI unitNo;
    tqgio((char*)"FILE ", &unitNo, &err);
    if (err) abort_prog(__LINE__, "tqgio", err);

    std::cout << "The thermochemical data will be read from the file "
              << "associated with unit "
              << unitNo
              << "\n\n";

    // Open ASCII thermochemical data-file for reading
    tqopna((char*)fileName.c_str(), unitNo, &err);
    if (err) abort_prog(__LINE__, "tqopna", err);

    // Read data-file
    tqrfil(&err);
    if (err) abort_prog(__LINE__, "tqrfil", err);

    // Close data-file
    tqclos(unitNo, &err);
    if (err) abort_prog(__LINE__, "tqclos", err);

    // Get default system units
    char dstr[TQSTRLEN];
    std::cout << "DEFAULT SYSTEM UNITS:\nQuantity    Unit\n";

    tqgsu((char*)"Pressure ", dstr, &err);
    if (err) abort_prog(__LINE__, "tqgsu", err);
    std::cout << "Pressure    " << dstr << "\n";

    tqgsu((char*)"Volume ", dstr, &err);
    if (err) abort_prog(__LINE__, "tqgsu", err);
    std::cout << "Volume      " << dstr << "\n";

    tqgsu((char*)"Temperature ", dstr, &err);
    if (err) abort_prog(__LINE__, "tqgsu", err);
    std::cout << "Temperature " << dstr << "\n";

    tqgsu((char*)"Energy ", dstr, &err);
    if (err) abort_prog(__LINE__, "tqgsu", err);
    std::cout << "Energy      " << dstr << "\n";

    tqgsu((char*)"Amount ", dstr, &err);
    if (err) abort_prog(__LINE__, "tqgsu", err);
    std::cout << "Amount      " << dstr << "\n\n";

    // Change system unit (default of 'Amount' is mol)
    // Here we change the unit for the quantity 'Amount' to gram, mainly so
    // that when we call 'tqstsc' firther down, we get the molecular mass
    // expressed in the unit g/mol
    tqcsu((char*)"Amount ", (char*)"gram ", &err);
    if (err) abort_prog(__LINE__, "tqcsu", err);
    std::cout << "Amount unit set to gram for masses below\n\n";

    // Get number of system components
    LI nSCom;
    tqnosc(&nSCom, &err);
    if (err) abort_prog(__LINE__, "tqnosc", err);

    std::cout << "SYSTEM COMPONENTS\nNumber of system compoents: "
              << nSCom
              << "\n";

    // Print the names of the system components and their masses
    std::cout << "No.  Name of component  Mass [g/mol]\n";
    for (LI i(1); i <= nSCom; ++i)
    {
        tqgnsc(i, dstr, &err);
        if (err) abort_prog(__LINE__, "tqgnsc", err);

        DB wMass;
        DB stoi[nSCom];
        tqstsc(i, stoi, &wMass, &err);
        if (err) abort_prog(__LINE__, "tqstsc", err);

        std::cout << i
                  << "    "
                  << std::setw(19)
                  << std::left
                  << dstr
                  << wMass
                  << "\n";
    }
    std::cout << "\n";

    // Get number of phases
    LI nPhase;
    char pName[TQSTRLEN];
    tqnop(&nPhase, &err);
    if (err) abort_prog(__LINE__, "tqnop", err);

    // Print the names of the phases and their model names
    std::cout << "PHASES\nNumber of phases: "
              << nPhase
              << "\nNo.  Name of phase  Model\n";

    for (LI i(1); i <= nPhase; ++i)
    {
        tqgnp(i, pName, &err);
        if (err) abort_prog(__LINE__, "tqgnp", err);
        tqmodl(i, dstr, &err);
        if (err) abort_prog(__LINE__, "tqmodl", err);

        std::cout << i
                  << "    "
                  << std::setw(15)
                  << std::left
                  << pName
                  << dstr
                  << "\n";
    }
    std::cout << "(PURE: Stoichiometric condensed phase, IDMX: Ideal mixing)"
              << "\n\n";

    // Get number of phase constituents of gas phase
    LI nPCon;
    tqnopc(1, &nPCon, &err);
    if (err) abort_prog(__LINE__, "tqnopc", err);

    // Print the names of the phase constituents and their Gibbs free energies
    // of formation
    std::cout << "PHASE CONSTITUENTS\nNumber of phase constituents of the "
              << "first gas phase: "
              << nPCon
              << "\nPh.  Name of constituent  Gibbs [J/mol]\n";

    char pCName[TQSTRLEN];
    for (LI i(1); i <= nPhase; ++i)
    {
        tqgnp(i, pName, &err);
        if (err) abort_prog(__LINE__, "tqgnp", err);
        tqnopc(i, &nPCon, &err);
        if (err) abort_prog(__LINE__, "tqnopc", err);

        for (LI j(1); j <= nPCon; ++j)
        {
            tqgnpc(i, j, pCName, &err);
            if (err) abort_prog(__LINE__, "tqgnpc", err);

            LI nValv;
            DB valv[25];
            tqgdat(i, j, (char*)"H ", 1, &nValv, valv, &err);
            if (err) abort_prog(__LINE__, "tqgdat", err);

            DB sEnthal(valv[0]);
            tqgdat(i, j, (char*)"S ", 1, &nValv, valv, &err);
            if (err) abort_prog(__LINE__, "tqgdat", err);
            DB sEntrop(valv[0]);

            // G = H - T*S
            DB sGibbs(sEnthal - (298.15 * sEntrop));
            std::cout << i
                      << "    "
                      << std::setw(19)
                      << std::left
                      << pCName
                      << "  "
                      << sGibbs
                      << "\n";
        }
    }
    std::cout << "\n";

    // Change system unit back to mol
    tqcsu((char*)"Amount ", (char*)"mol ", &err);
    if (err) abort_prog(__LINE__, "tqcsu", err);
    std::cout << "Amount unit set back to mol" << "\n\n";

    // Make sure everything can be used as an incoming species
    std::cout << "Mixture phase  Phase constituent  Ok/Not permitted\n";

    LI nPerm(0);
    for (LI i(1); i <= nPhase; ++i)
    {
        tqgnp(i, pName, &err);
        if (err) abort_prog(__LINE__, "tqgnp", err);
        tqnopc(i, &nPCon, &err);
        if (err) abort_prog(__LINE__, "tqnopc", err);

        // We only need to check mixture phases
        if (nPCon > 1)
        {
            for (LI j(1); j <= nPCon; ++j)
            {
                tqgnpc(i, j, pCName, &err);
                if (err) abort_prog(__LINE__, "tqgnpc", err);

                // Check whether the current constituent is permitted to be
                // used as incoming species
                LI inPCIS;
                tqpcis(i, j, &inPCIS, &err);
                if (err) abort_prog(__LINE__, "tqpcis", err);

                // Print a table entry, and if it is not permitted, increase
                // the total number of such phase constituents by 1
                std::cout << std::setw(15)
                          << std::left
                          << pName
                          << std::setw(19)
                          << std::left
                          << pCName;
                if (inPCIS == 0)
                {
                    ++nPerm;
                    std::cout << "Not peritted\n";
                } else {
                    std::cout << "Ok\n";
                }
            }
        }
    }

    // If there were phase constituents which cannot be used as incoming
    // species, print a note
    if (nPerm > 0)
        std::cout << "Note: "
                  << nPerm
                  << " phase constituent(s) is/are not "
                  << "permitted as incoming species\n";
}


void open_file(std::string logDir, LI& err)
{
    // Truncate file containing molar masses
    std::ofstream molMassFile(logDir + "/mole_masses.dat",
                              std::ofstream::trunc);

    // Redirect the results (output of ChemApp routines like 'tqcel' and
    // 'tqshow') to an output file for viewing
    tqcio((char*)"LIST ", 3, &err);
    if (err) abort_prog(__LINE__, "tqcio", err);

    std::string resultsFileName(logDir + "/results.dat");
    tqopen((char*)resultsFileName.c_str(), 3, &err);
    if (err) abort_prog(__LINE__, "tqopen", err);
}


std::vector<double> amounts(std::vector<double> concs, std::string waterRole,
                            DB temperature, DB pressure,
                            std::string logDir, LI& err)
{
    // Set temperature and pressure
    LI numCon;
    tqsetc((char*)"T ", 0, 0, temperature, &numCon, &err);
    if (err) abort_prog(__LINE__, "tqsetc", err);
    tqsetc((char*)"P ", 0, 0, pressure, &numCon, &err);
    if (err) abort_prog(__LINE__, "tqsetc", err);

    // Get number of phases
    LI nPhase;
    tqnop(&nPhase, &err);
    if (err) abort_prog(__LINE__, "tqnop", err);

    // Set initial concentrations in ChemApp
    int indexConcs(0);
    LI nPCon;
    for (LI i(1); i <= nPhase; ++i)
    {
        // Find number of constituents in phase
        tqnopc(i, &nPCon, &err);
        if (err) abort_prog(__LINE__, "tqnopc", err);

        for (LI j(1); j <= nPCon; ++j)
        {
            tqsetc((char*)"ia ", i, j, static_cast<DB>(concs[indexConcs]),
                   &numCon, &err);
            if (err) abort_prog(__LINE__, "tqsetc", err);

            ++indexConcs;
        }
    }

    // Use 'tqwstr' to write header to the output file
    tqwstr((char*)"LIST ", (char*)"RESULTS ", &err);
    if (err) abort_prog(__LINE__, "tqwstr", err);

    // Calculate and write equilibrium to the output file
    // The "  " option means there will be no target variable estimate used
    // to begin the simulator at, meaning 'estim" istn't used
    DB estim[1] = { 0 };
    tqcel((char*)" ", 0, 0, estim, &err);
    if (err) abort_prog(__LINE__, "tqcel", err);

    // Determine new concentration of water for renormalization of
    // concentrations of constituents.
    // If water takes part in the reaction, collect corresponding indices
    std::vector< std::pair<LI, LI> > indicesWater;
    for (LI i(1); i <= nPhase; ++i)
    {
        // Find number of constituents in phase
        tqnopc(i, &nPCon, &err);
        if (err) abort_prog(__LINE__, "tqnopc", err);

        for (LI j(1); j <= nPCon; ++j)
        {
            // Determine name of constituent
            char pCName[TQSTRLEN];
            tqgnpc(i, j, pCName, &err);
            if (err) abort_prog(__LINE__, "tqgnpc", err);

            // Check if constituent is water and if collect corresponding
            // indices
            char h2oStr[TQSTRLEN] = "H2O";
            if (!std::strcmp(pCName, h2oStr))
                indicesWater.push_back(std::make_pair(i, j));
        }
    }

    // If part of reaction, get concentration of water
    DB waterConc(1.0);
    if (waterRole == "reactant")
        waterConc = 0.0;
    // Throw exception, if role of water is given as "solvent", but water
    // was given in reaction
    if (waterRole == "solvent" && indicesWater.size() != 0)
    {
        std::string errorStr("Role of water was given inconsistent to " \
                             "appearence of water as reagent in reaction");
        throw ErrorStr(errorStr);
    }
    if (indicesWater.size() != 0)
    {
        for (std::pair<LI, LI> indices:indicesWater)
        {
            DB phaseWaterConc;
            tqgetr((char*)"A ", indices.first, indices.second,
                   &phaseWaterConc, &err);
            if (err) abort_prog(__LINE__, "tqgetr", err);

            if (waterRole == "reactant" || waterRole == "product")
                waterConc += phaseWaterConc;
            else
            {
                std::string errorStr("role of water was given " \
                                     "inconsistent to appearence of " \
                                     "water as reagent in reaction");
                throw ErrorStr(errorStr);
            }
        }
    }

    // Collect amounts to save them in 'concs'
    indexConcs = 0;
    for (LI i(1); i <= nPhase; ++i)
    {
        tqnopc(i, &nPCon, &err);
        if (err) abort_prog(__LINE__, "tqnopc", err);

        for (LI j(1); j <= nPCon; ++j)
        {
            // Get amounts of constituents
            DB constiAmount;
            tqgetr((char*)"A ", i, j, &constiAmount, &err);
            if (err) abort_prog(__LINE__, "tqgetr", err);

            // Normalize constituent concentration to new found water
            // concentration and save to 'concs' for return
            concs[indexConcs] = static_cast<double>(constiAmount / waterConc);
            ++indexConcs;
        }
    }

    return concs;
}


void close_file(LI& err)
{
    // Close file containing redirected results of ChemApp routines like
    // 'tqcel' and 'tqshow'
    tqclos(3, &err);
    if (err) abort_prog(__LINE__, "tqclos", err);

    // Redirect the output back to unit number 6 (standard output unit)
    tqcio((char*)"LIST ", 6, &err);
    if (err) abort_prog(__LINE__, "tqcio", err);
}


PYBIND11_MODULE(ChemApp, m)
{
    m.doc() = "Calculate equlibrium amounts of nucleobases using a C++ port "
              "of FORTRAN library ChemApp.";

    m.def("start", &start, "Start ChemApp.",
          "err"_a,
          py::call_guard<py::scoped_ostream_redirect,
                         py::scoped_estream_redirect>());
    m.def("read_data", &read_data, "Read input file into ChemApp.",
          "logDir"_a, "err"_a,
          py::call_guard<py::scoped_ostream_redirect,
                         py::scoped_estream_redirect>());
    m.def("open_file", &open_file, "Open log files.",
          "logDir"_a, "err"_a,
          py::call_guard<py::scoped_ostream_redirect,
                         py::scoped_estream_redirect>());
    m.def("amounts", &amounts, "Calculate equilibrium amounts of constituents "
          "and return as array of floats.",
          "waterRole"_a, "concs"_a, "temperature"_a, "pressure"_a, "logDir"_a,
          "err"_a,
          py::call_guard<py::scoped_ostream_redirect,
                         py::scoped_estream_redirect>());
    m.def("close_file", &close_file, "Close log files.",
          "err"_a,
          py::call_guard<py::scoped_ostream_redirect,
                         py::scoped_estream_redirect>());
}
