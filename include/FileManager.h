#ifndef FILEMANAGER_H
#define FILEMANAGER_H

#include "TFile.h"
#include "Notifier.h"

class FileManager
{
public:
    FileManager();
    ~FileManager(void);
    void ProcessInputFiles(int argc, char **argv);
    void OpenRootFile(std::string file_name);
    void LoadInternalCalibration(TFile *file);
    void AutoFileDetect(std::string file_name);

    TChain *GetChain()
    {
        return gChain;
    }

    TChain *gChain = NULL;

private:
    Notifier *notifier;
};

#endif