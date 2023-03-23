#ifndef SORT_ANALYSIS_TREES_H
#define SORT_ANALYSIS_TREES_H

#include "Notifier.h"

int main(int argc, char **argv);
void AutoFileDetect(std::string file_name);
void OpenRootFile(std::string file_name);
void LoadInternalCalibration();
void PrintUsage(char *argv[]);
void InitGRSISort();

std::string grsi_path;
Notifier *notifier = new Notifier;

#endif
