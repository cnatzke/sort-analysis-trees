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
#include "TFile.h" // needed for GetRunNumber
#include "TGRSIUtilities.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "HistogramManager.h"
#include "SortAnalysisTrees.h"

int main(int argc, char **argv)
{
    std::string compton_limits_filepath = "/data1/S9038/current-sort/data/histograms/compton-algorithm/compton_limits.csv";

    if (argc == 1) { // no inputs given
        PrintUsage(argv);
        return 0;
    }

    InitGRSISort();

    for (auto i = 1; i < argc; i++) AutoFileDetect(argv[i]);

    if (!gChain) std::cout << "No gChain found" << std::endl;
    if (!gChain->GetEntries()) std::cout << "Found gChain, but no entries retrieved" << std::endl;

    if (!gChain || !gChain->GetEntries()) {
        std::cerr << "Failed to find anything. Exiting" << std::endl;
        return 1;
    }

    // Specifiy crosstalk corrections
    TGRSIOptions::AnalysisOptions()->SetCorrectCrossTalk(false);

    HistogramManager * hist_man = new HistogramManager(compton_limits_filepath);
    // turn on multiplicity filter
    //hist_man->SetMultiplicityFilter(true);

    // fill histograms
    hist_man->MakeHistograms(gChain);
    hist_man->WriteHistogramsToFile();

    delete hist_man;
    std::cout << "Exiting ... " << std::endl;
    return 0;
} // main()

/************************************************************//**
 * Opens Root files
 *
 * @param file_name Analysis file path
 ***************************************************************/
void OpenRootFile(std::string file_name)
{
    TFile f(file_name.c_str());
    if (f.Get("AnalysisTree")) {
        if (!gChain) {
            gChain = new TChain("AnalysisTree");
            notifier->AddChain(gChain);
            gChain->SetNotify(notifier);
        }
        gChain->Add(file_name.c_str());
        //std::cout << "Added: " << file_name << std::endl;
    }
} // end OpenRootFile

/******************************************************************************
 * Determines input file type
 *
 * @param file_name  Input file
 *****************************************************************************/
void AutoFileDetect(std::string file_name)
{
    size_t dot_pos = file_name.find_last_of('.');
    std::string ext = file_name.substr(dot_pos + 1);

    if (ext == "root") {
        OpenRootFile(file_name);
    }
    else if (ext == "cal") {
        notifier->AddCalFile(file_name);
    } else {
        std::cerr << "Discarding unknown file: " << file_name.c_str() << std::endl;
    }
} // AutoFileDetect()

/******************************************************************************
 * Initializes GRSISort libraries and stuff
 *
 *****************************************************************************/
void InitGRSISort(){
    // makes time retrival happy and loads GRSIEnv
    grsi_path = getenv("GRSISYS");
    if(grsi_path.length() > 0) {
        grsi_path += "/";
    }
    grsi_path += ".grsirc";
    gEnv->ReadFile(grsi_path.c_str(), kEnvChange);

    TParserLibrary::Get()->Load();
} // end InitGRSISort

/******************************************************************************
 * Prints usage message and version
 *****************************************************************************/
void PrintUsage(char* argv[]){
    std::cerr << argv[0] << " Version: " << SortAnalysisTrees_VERSION_MAJOR << "." << SortAnalysisTrees_VERSION_MINOR << "\n"
              << "usage: " << argv[0] << " calibration_file analysis_tree [analysis_tree_2 ... ] calibration_file\n"
              << " calibration_file:       calibration file (must end with .cal)\n"
              << " analysis_tree:          analysis tree to process (must end with .root)"
              << std::endl;
} // end PrintUsage
