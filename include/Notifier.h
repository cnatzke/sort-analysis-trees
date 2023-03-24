#ifndef NOTIFIER_H
#define NOTIFIER_H

#include "TChain.h"
#include "TPPG.h"

class Notifier : public TObject
{
public:
	Notifier()
	{
	}
	~Notifier()
	{
	}

	void AddChain(TChain *chain)
	{
		fChain = chain;
	}
	void AddRootFile(std::string name)
	{
		RootFiles.push_back(name);
	}
	void AddInfoFile(std::string name)
	{
		InfoFiles.push_back(name);
	}
	void AddCalFile(std::string name)
	{
		CalFiles.push_back(name);
	}
	void InternalCalFile(bool flag)
	{
		internal_cal_file_flag = flag;
	}
	bool GetInternalCalFlag()
	{
		return internal_cal_file_flag;
	}

	bool Notify()
	{
		printf("%s loaded.\n", fChain->GetCurrentFile()->GetName());
		ppg = (TPPG *)fChain->GetCurrentFile()->Get("TPPG");

		if (CalFiles.size() > 0)
		{
			TChannel::ReadCalFile(CalFiles.at(0).c_str());
		}
		else if (internal_cal_file_flag)
		{
			std::cout << "Loaded calibration from ROOT file" << std::endl;
		}
		else
		{
			std::cout << "No calibration file loaded." << std::endl;
		}

		return true;
	}

	// TGRSIRunInfo *info = NULL;
	TPPG *ppg = NULL;
	TList *outList = NULL;

	bool internal_cal_file_flag = true;

private:
	TChain *fChain;
	std::vector<std::string> RootFiles;
	std::vector<std::string> CalFiles;
	std::vector<std::string> InfoFiles;

}; // End Notifier class
#endif
