//////////////////////////////////////////////////////////////////////////////////
// Recovers intra-clover Compton Scatters
//
// Author: Connor Natzke (cnatzke@triumf.ca)
// Creation Date: 2022-03-15
// Last Update:   2022-03-15
//////////////////////////////////////////////////////////////////////////////////
#include "ComptonRecovery.h"

#include <fstream>
#include "csv.h"
#include "TMath.h"

/************************************************************//**
 * Constructor
 ***************************************************************/
ComptonRecovery::ComptonRecovery(std::string file_path)
{
    _energy_resolution = 0.0023; // 0.23%
    // read in Compton angular limits
    ReadInComptonLimits(file_path);
} // end Constructor

/************************************************************//**
 * Destructor
 ***************************************************************/
ComptonRecovery::~ComptonRecovery(void)
{
} // end Destructor

/************************************************************//**
 * Checks if the angle and energy between two photons can
 * be a Compton scatter
 *
 * @param angle Angle between two photons (radians)
 * @param energy_1 Energy of first photon (keV)
 * @param energy_2 Energy of second photon (keV)
 *****************************************************************************/
bool ComptonRecovery::ComptonScatterCandidate(int angle_index, float energy_1, float energy_2)
{
    int verbose = 0;
    float energy_total = energy_1 + energy_2;
    float compton_angle_ftb = _compton_limit_vec_ftb.at(angle_index);
    float compton_angle_btf = _compton_limit_vec_btf.at(angle_index);

    double scattered_photon_energy_upper = ComptonScatterFunction(compton_angle_ftb, energy_total * (1.0 + 2 * _energy_resolution));
    double scattered_photon_energy_lower = ComptonScatterFunction(compton_angle_btf, energy_total * (1.0 - 2 * _energy_resolution));
    double second_photon_energy_upper = energy_total - scattered_photon_energy_lower;
    double second_photon_energy_lower = energy_total - scattered_photon_energy_upper;

    if (((scattered_photon_energy_lower < energy_1) && (energy_1 < scattered_photon_energy_upper)) && ((second_photon_energy_lower < energy_2) && (energy_2 < second_photon_energy_upper))) {
        if (verbose > 0) {
            std::cout << "[" << scattered_photon_energy_lower << ", " << scattered_photon_energy_upper << "] " << energy_1 << std::endl;
            std::cout << "[" << second_photon_energy_lower << ", " << second_photon_energy_upper << "] " << energy_2 << std::endl;
            std::cout << "-----------" << std::endl;
        }
        return true;
    }
    else if (((scattered_photon_energy_lower < energy_2) && (energy_2 < scattered_photon_energy_upper)) && ((second_photon_energy_lower < energy_1) && (energy_1 < second_photon_energy_upper))) {
        if (verbose > 0) {
            std::cout << "[" << scattered_photon_energy_lower << ", " << scattered_photon_energy_upper << "] " << energy_2 << std::endl;
            std::cout << "[" << second_photon_energy_lower << ", " << second_photon_energy_upper << "] " << energy_1 << std::endl;
            std::cout << "-----------" << std::endl;
        }
        return true;
    } else {
        return false;
    }
} // end ComptonScatterCandidate()

/************************************************************//**
 * Recovers Compton scatter event
 *
 * @param hit1 First GRIFFIN hit
 * @param hit2 Second GRIFFIN hit
 *****************************************************************************/
void ComptonRecovery::RecoverComptonScatter(TGriffinHit * hit1, TGriffinHit * hit2)
{
} // end RecoverComptonScatter()

/************************************************************//**
 * Returns the position of highest energy deposit
 *
 * @param hit1 First GRIFFIN hit
 * @param hit2 Second GRIFFIN hit
 *****************************************************************************/
bool ComptonRecovery::FirstHitHigh(TGriffinHit * hit1, TGriffinHit * hit2)
{
    // assign reconstucted location to highest energy deposit
    if (hit1->GetEnergy() > hit2->GetEnergy()) {
        return true;
    } else {
        return false;
    }
} // end FirstHitHigh()

/************************************************************//**
 * Returns the time of the first hit
 *
 * @param hit1 First GRIFFIN hit
 * @param hit2 Second GRIFFIN hit
 *****************************************************************************/
double ComptonRecovery::GetReconstructedTime(TGriffinHit * hit1, TGriffinHit * hit2)
{
    // assign time as time of the first hit
    double delta_t = (hit1->GetTime() - hit2->GetTime());

    if (delta_t > 0){
        return hit2->GetTime();
    } else {
        return hit1->GetTime();
    }
} // end GetReconstructedTime

/************************************************************//**
 * Compton scatter formula
 *
 * @param angle Compton scatter angle (rad)
 * @param energy Energy of initial photon (keV)
 *****************************************************************************/
double ComptonRecovery::ComptonScatterFunction(double angle, float energy)
{
    float electron_rest_mass_energy = 511.; //keV
    return energy / (1 + (energy / electron_rest_mass_energy) * (1 - TMath::Cos(angle)));
} // end ComptonScatter()

/************************************************************//**
 * Reads in csv file with Compton scattering bands
 *
 * @param filepath Filepath to csv file
 ***************************************************************/
void ComptonRecovery::ReadInComptonLimits(std::string filepath)
{
    std::fstream _compton_limits_file;

    _compton_limits_file.open(filepath, std::ios_base::in);
    // if bg file doesn't exist, throw error
    if (!_compton_limits_file) {
        _compton_limits_file.close();
        std::cout << "Could not open " << filepath << ", exiting." << std::endl;
        exit(EXIT_FAILURE);
    } else {
        std::cout << "Found accepted Compton scatter angles file: " << filepath << "\n" << std::endl;
        io::CSVReader<5> in(filepath);
        in.read_header(io::ignore_extra_column, "angular_bin", "compton_limit_ftb", "compton_limit_btf", "angle_diff_ftb", "angle_diff_btf");
        float angular_bin; float compton_limit_ftb; float compton_limit_btf; float angle_diff_ftb; float angle_diff_btf;
        while(in.read_row(angular_bin, compton_limit_ftb, compton_limit_btf, angle_diff_ftb, angle_diff_btf)) {
            _angular_bin_vec.push_back(angular_bin);
            _compton_limit_vec_ftb.push_back(compton_limit_ftb);
            _compton_limit_vec_btf.push_back(compton_limit_btf);
        }
        _compton_limits_file.close();
        // remove first bin of filled vectors sinces it's the zero angle
        _angular_bin_vec.erase(_angular_bin_vec.begin());
        _compton_limit_vec_ftb.erase(_compton_limit_vec_ftb.begin());
        _compton_limit_vec_btf.erase(_compton_limit_vec_btf.begin());
    }

} // end ReadInComptonLimits
