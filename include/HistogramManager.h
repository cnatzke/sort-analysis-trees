#ifndef HISTOGRAM_MANAGER_H
#define HISTOGRAM_MANAGER_H

#include <map>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TGriffin.h"
#include "TGriffinBgo.h"
#include "TChain.h"
#include "TVector3.h"
#include "TMath.h"
#include "ComptonRecovery.h"

class HistogramManager
{
std::string _compton_limits_file;
bool _multiplicity_filter = false;
int _multiplicity_limit = 2;
int _event_number;
int _prompt_time_max;
int _bg_time_min;
float _detector_radius;

// Histogram vectors
std::map<std::string, TH1D*> hist_1D;
std::map<std::string, TH2D*> hist_2D;

// energy vectors
std::vector<float> crystal_energy_vec;
std::vector<float> addback_energy_vec;
//std::vector<float> unsup_crystal_energy_vec;
//std::vector<float> reconstructed_crystal_energy_vec;
std::vector<float> reconstructed_addback_energy_vec;
// position vectors
std::vector<TVector3> crystal_pos_vec;
std::vector<TVector3> addback_pos_vec;
//std::vector<TVector3> unsup_crystal_pos_vec;
//std::vector<TVector3> reconstructed_crystal_pos_vec;
std::vector<TVector3> reconstructed_addback_pos_vec;
// time vectors
std::vector<double> crystal_time_vec;
std::vector<double> addback_time_vec;
//std::vector<double> unsup_crystal_time_vec;
//std::vector<double> reconstructed_crystal_time_vec;
std::vector<double> reconstructed_addback_time_vec;
// detector id vectors
std::vector<int> crystal_id_vec;
std::vector<int> addback_id_vec;
std::vector<int> crystal_clover_id_vec;
//std::vector<int> unsup_crystal_id_vec;
//std::vector<int> reconstructed_crystal_id_vec;
std::vector<int> reconstructed_addback_id_vec;
// misc vectors
std::vector<float> crystal_kvalue_vec;

public:
HistogramManager(std::string compton_limits_file);
~HistogramManager(void);
void MakeHistograms(TChain *inputChain);
void FillHistograms(TChain *gChain);
void InitializeHistograms(int verbose = 0);
void WriteHistogramsToFile();
void SetMultiplicityFilter(bool val) {
    _multiplicity_filter = val;
};
bool GetMultiplicityFilter() {
    return _multiplicity_filter;
};

private:
void PreProcessData(ComptonRecovery * comp_check);
int GetAngleIndex(double angle, std::vector<double> vec);
int GetClosest(int val1, int val2, std::vector<double> vec, double target);
bool IsInSlice(double delta_t, int prompt_time);

double degree_to_rad = TMath::Pi() / 180.;
double rad_to_degree = 180. / TMath::Pi();
TGriffin *fGrif = NULL;
TGriffinBgo *fGriffinBgo = NULL;

int num_crystals = 64;
std::vector<double> angle_combinations_vec = {15.442, 21.9054, 29.1432, 33.1433, 38.382, 44.57, 47.4453, 48.7411, 51.4734, 55.1704, 59.9782, 60.1024, 62.3396, 62.4924, 63.4231, 68.9567, 71.4314, 73.3582, 73.6291, 75.7736, 80.9423, 81.5464, 83.8936, 86.868, 88.9658, 91.0342, 93.132, 96.1064, 98.4536, 99.0577, 104.226, 106.371, 106.642, 108.569, 111.043, 116.577, 117.508, 117.66, 119.898, 120.022, 124.83, 128.527, 131.259, 132.555, 135.43, 141.618, 146.857, 150.857, 158.095, 164.558, 180.0};

double offsets[64];
double gains[64];


};

#endif
