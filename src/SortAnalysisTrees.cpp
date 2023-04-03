//////////////////////////////////////////////////////////////////////////////////
// Creates histograms from unpacked analysis trees
//
// Author:        Connor Natzke (cnatzke@triumf.ca)
// Creation Date: June 2020
// Last Update:   June 2020
// Usage:
//
//////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "TGRSIUtilities.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "FileManager.h"
#include "HistogramManager.h"
#include "SortAnalysisTrees.h"

int main(int argc, char **argv)
{
    bool use_compton_alg = false;

    if (argc == 1)
    { // no inputs given
        PrintUsage(argv);
        return 0;
    }

    InitGRSISort();

    FileManager *file_manager = new FileManager();
    file_manager->ProcessInputFiles(argc, argv);
    event_chain = file_manager->GetChain();

    if (!event_chain)
        std::cout << "No event_chain found" << std::endl;
    if (!event_chain->GetEntries())
        std::cout << "Found event_chain, but no entries retrieved" << std::endl;

    if (!event_chain || !event_chain->GetEntries())
    {
        std::cerr << "Failed to find anything. Exiting" << std::endl;
        return 1;
    }

    // Specifiy crosstalk corrections
    TGRSIOptions::AnalysisOptions()->SetCorrectCrossTalk(false);

    HistogramManager *hist_man;
    if (use_compton_alg)
    {
        std::string compton_limits_filepath = "./compton_limits_145mm.csv";
        hist_man = new HistogramManager(compton_limits_filepath);
    }
    else
    {
        hist_man = new HistogramManager();
    }
    // turn on multiplicity filter
    // hist_man->SetMultiplicityFilter(true);

    // fill histograms
    hist_man->MakeHistograms(event_chain);
    hist_man->WriteHistogramsToFile();

    delete file_manager;
    delete hist_man;
    std::cout << "Exiting ... " << std::endl;

    return 0;
} // main()

/******************************************************************************
 * Initializes GRSISort libraries and stuff
 *
 *****************************************************************************/
void InitGRSISort()
{
    // makes time retrival happy and loads GRSIEnv
    grsi_path = getenv("GRSISYS");
    if (grsi_path.length() > 0)
    {
        grsi_path += "/";
    }
    grsi_path += ".grsirc";
    gEnv->ReadFile(grsi_path.c_str(), kEnvChange);

    TParserLibrary::Get()->Load();
} // end InitGRSISort

/******************************************************************************
 * Prints usage message and version
 *****************************************************************************/
void PrintUsage(char *argv[])
{
    std::cerr << argv[0] << "\n"
              << "usage: " << argv[0] << " analysis_tree [analysis_tree_2 ... ] calibration_file\n"
              << " analysis_tree:          analysis tree to process (must end with .root)\n"
              << " calibration_file:       calibration file (optional, must end with .cal)"
              << std::endl;
} // end PrintUsage
