//////////////////////////////////////////////////////////////////////////////////
// Manages file IO
//
// Author:        Connor Natzke (cnatzke@triumf.ca)
// Creation Date: Mar 2023
// Last Update:   Mar 2023
// Usage:
//
//////////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "TParserLibrary.h"
#include "FileManager.h"

/**************************************************************
 * Constructors
 ***************************************************************/
FileManager::FileManager()
{
    notifier = new Notifier();
} // end Constructor

/****************************************************************
 * Destructor
 ***************************************************************/
FileManager::~FileManager(void)
{
} // end Destructor

/******************************************************************************
 * Opens input files and determines what calibration to use
 *
 * @param file_name  Input file
 *****************************************************************************/
void FileManager::ProcessInputFiles(int argc, char **argv)
{

    for (auto i = 1; i < argc; i++)
    {
        AutoFileDetect(argv[i]);
    }

} // ProcessInputFiles

/******************************************************************************
 * Determines input file type
 *
 * @param file_name  Input file
 *****************************************************************************/
void FileManager::AutoFileDetect(std::string file_name)
{
    size_t dot_pos = file_name.find_last_of('.');
    std::string ext = file_name.substr(dot_pos + 1);

    if (ext == "root")
    {
        OpenRootFile(file_name);
    }
    else if (ext == "cal")
    {
        notifier->InternalCalFile(false);
        notifier->AddCalFile(file_name);
    }
    else
    {
        std::cerr << "Discarding unknown file: " << file_name.c_str() << std::endl;
    }
} // AutoFileDetect()

/**************************************************************
 * Opens Root files
 *
 * @param file_name Analysis file path
 ***************************************************************/
void FileManager::OpenRootFile(std::string file_name)
{
    TFile f(file_name.c_str());
    if (f.Get("AnalysisTree"))
    {
        if (!gChain)
        {
            gChain = new TChain("AnalysisTree");
            notifier->AddChain(gChain);
            gChain->SetNotify(notifier);
        }
        gChain->Add(file_name.c_str());
        // std::cout << "Added: " << file_name << std::endl;
    }
    LoadInternalCalibration(&f);
} // end OpenRootFile

/**************************************************************
 * Reads calibration parameters from the root file
 *
 ***************************************************************/
void FileManager::LoadInternalCalibration(TFile *file)
{
    TChannel::ReadCalFromFile(file);
}