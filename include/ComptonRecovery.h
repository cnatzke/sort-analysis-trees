#ifndef COMPTON_RECOVERY_H
#define COMPTON_RECOVERY_H

#include <vector>
#include <iostream>
#include "TGriffin.h"

class ComptonRecovery
{
float _energy_resolution;
std::vector<float> _angular_bin_vec;
std::vector<float> _compton_limit_vec_ftb;
std::vector<float> _compton_limit_vec_btf;

public:
    ComptonRecovery(std::string file_path);
    ~ComptonRecovery(void);
    bool ComptonScatterCandidate(int angle_index, float energy_1, float energy_2);
    void RecoverComptonScatter(TGriffinHit * hit1, TGriffinHit * hit2);
    bool FirstHitHigh(TGriffinHit * hit1, TGriffinHit * hit2);
    double GetReconstructedTime(TGriffinHit * hit1, TGriffinHit * hit2);
    double ComptonScatterFunction(double angle, float energy);
    void ReadInComptonLimits(std::string filepath);

};

#endif
