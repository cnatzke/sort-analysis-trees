//////////////////////////////////////////////////////////////////////////////////
// Creates and fills histograms
//
// Author:          Connor Natzke (cnatzke@triumf.ca)
// Creation Date: Wednesday July 29, 2020	T15:22:33-07:00
// Last Update:   Wednesday July 29, 2020	T15:22:45-07:00
// Usage:
//
//////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <map>
#include "csv.h"
#include "HistogramManager.h"
#include "progress_bar.h"
#include "LoadingMessenger.h"


/************************************************************//**
 * Creates and Fills histograms
 *
 * @param verbose Verbosity level
 ***************************************************************/
void HistogramManager::MakeHistograms(TChain *input_chain)
{
    int verbosity = 1;
    std::string compton_limits_filepath = "/data1/S9038/current-sort/data/histograms/compton-algorithm/compton_limits.csv";
    // reads in Compton limits from file
    ReadInComptonLimits(compton_limits_filepath);
    // create histograms
    InitializeHistograms(verbosity);
    FillHistograms(input_chain);

} // GenerateHistogramFile()


/************************************************************//**
 * Reads in csv file with Compton scattering bands
 *
 * @param filepath Filepath to csv file
 ***************************************************************/
void HistogramManager::ReadInComptonLimits(std::string filepath)
{
    std::fstream compton_limits_file;

    compton_limits_file.open(filepath, std::ios_base::in);
    // if bg file doesn't exist, throw error
    if (!compton_limits_file) {
        compton_limits_file.close();
        std::cout << "Could not open " << filepath << ", exiting." << std::endl;
        exit(EXIT_FAILURE);
    } else {
        std::cout << "Found accepted Compton scatter angles file: " << filepath << std::endl;
        io::CSVReader<5> in(filepath);
        in.read_header(io::ignore_extra_column, "angular_bin", "compton_limit_ftb", "compton_limit_btf", "angle_diff_ftb", "angle_diff_btf");
        float angular_bin; float compton_limit_ftb; float compton_limit_btf; float angle_diff_ftb; float angle_diff_btf;
        while(in.read_row(angular_bin, compton_limit_ftb, compton_limit_btf, angle_diff_ftb, angle_diff_btf)) {
            angular_bin_vec.push_back(angular_bin);
            compton_limit_vec_ftb.push_back(compton_limit_ftb);
            compton_limit_vec_btf.push_back(compton_limit_btf);
        }
        compton_limits_file.close();
        // remove first bin of filled vectors sinces it's the zero angle
        angular_bin_vec.erase(angular_bin_vec.begin());
        compton_limit_vec_ftb.erase(compton_limit_vec_ftb.begin());
        compton_limit_vec_btf.erase(compton_limit_vec_btf.begin());
    }

} // end ReadInComptonLimits

/************************************************************//**
 * Initializes histograms to be filled
 *
 * @param verbose Verbosity level
 ***************************************************************/
void HistogramManager::InitializeHistograms(int verbose)
{
    int energy_bins_max = 3000;
    int energy_bins_min = 0;

    if (verbose > 0 ) std::cout << "Creating 1D histograms ... " << std::endl;

    // 1D Histograms
    hist_1D["delta_t"] = new TH1D("delta_time", "#Delta t", energy_bins_max, energy_bins_min, energy_bins_max);
    hist_1D["k_value"] = new TH1D("k_value", "#Delta t", 1000, 0, 1000);
    hist_1D["gamma_energy"] = new TH1D("gamma_energy", "gamma singles", energy_bins_max, energy_bins_min, energy_bins_max);

    // 2D Histograms
    if (verbose > 0 ) std::cout << "Creating 2D histograms ... " << std::endl;
    hist_2D["singles_channel"] = new TH2D("singles_channel", "#gamma singles vs channel", 70, 0, 70, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["gg_singles"] = new TH2D("gg_singles", "", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["sum_energy_angle"] = new TH2D("sum_energy_angle", "Sum Energy vs angle index", 70, 0, 70, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["sum_energy_angle_tr"] = new TH2D("sum_energy_angle_tr", "Sum Energy vs angle index (time random);angle index;sum energy [keV]", 70, 0, 70, energy_bins_max, energy_bins_min, energy_bins_max);

    // individual histograms for each angular bin
    for (unsigned int i = 0; i < angle_combinations_vec.size(); i++) {
        hist_2D_prompt[Form("index_%02i_sum", i)] = new TH2D(Form("index_%02i_sum", i), ";sum energy [keV];#gamma_1 energy [keV]", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
        // Time randoms
        hist_2D_tr[Form("index_%02i_sum_tr", i)] = new TH2D(Form("index_%02i_sum_tr", i), ";sum energy [keV];#Deltat [ns]", energy_bins_max, energy_bins_min, energy_bins_max, 600, 400, 1000);
        hist_2D_tr[Form("index_%02i_sum_tr_avg", i)] = new TH2D(Form("index_%02i_sum_tr_avg", i), ";sum energy [keV];#gamma_1 energy [keV]", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
        // Compton algorithm histograms
        hist_2D_comp_alg_acc[Form("index_%02i", i)] = new TH2D(Form("index_%02i_comp_acc", i), Form("Compton Algorithm Accepted %.2f;sum energy [keV];#gamma_1 energy [keV]", angle_combinations_vec.at(i)), energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
        hist_2D_comp_alg_rej[Form("index_%02i", i)] = new TH2D(Form("index_%02i_comp_rej", i), Form("Compton Algorithm Rejected %.2f;sum energy [keV];#gamma_1 energy [keV]", angle_combinations_vec.at(i)), energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
    }

} // InitializeHistograms()

/************************************************************//**
 * Fills histograms
 *
 * @param gChain Data chain
 ***************************************************************/
void HistogramManager::FillHistograms(TChain *gChain)
{

    int prompt_time_max = 30;     // ns
    int bg_time_min = 500;     // ns

    if (gChain->FindBranch("TGriffin")) {
        gChain->SetBranchAddress("TGriffin", &fGrif);
        if (fGrif != NULL) {
            std::cout << "Successfully found TGriffin branch" << std::endl;
        } else {
            std::cout << "Could not find TGriffin branch ... exiting" << std::endl;
        }
    }
    if (gChain->FindBranch("TGriffinBgo")) {
        gChain->SetBranchAddress("TGriffinBgo", &fGriffinBgo);
        if (fGriffinBgo != NULL) {
            std::cout << "Successfully found TGriffinBgo branch" << std::endl;
        } else {
            std::cout << "Could not find TGriffinBgo branch ... exiting" << std::endl;
        }
    }

    //DisplayLoadingMessage();
    LoadingMessenger load_man;
    load_man.DisplayLoadingMessage();

    long analysis_entries = gChain->GetEntries();
    ProgressBar progress_bar(analysis_entries, 70, '=', ' ');
    for (auto i = 0; i < analysis_entries; i++) {
        gChain->GetEntry(i);

        // Applies secondary energy calculation
        PreProcessData();

        // Filling histograms
        if (energy_vec.size() > 0) {
            for (unsigned int g1 = 0; g1 < energy_vec.size(); ++g1) {
                hist_1D["k_value"]->Fill(kvalue_vec.at(g1));
                hist_1D["gamma_energy"]->Fill(energy_vec.at(g1));
                hist_2D["singles_channel"]->Fill(detector_vec.at(g1), energy_vec.at(g1));
                // for(unsigned int g2 = 0; g2 < energy_vec.size(); ++g2) { // Makes matrices symmetric
                for(unsigned int g2 = g1 + 1; g2 < energy_vec.size(); ++g2) {                 // Makes matrices assymmetric
                    if (g1 == g2) continue;

                    double angle_rad = pos_vec.at(g1).Angle(pos_vec.at(g2));
                    double angle = angle_rad * rad_to_degree;
                    if (angle < 0.0001 || angle > 180.0) {
                        continue;
                    }

                    int angle_index = GetAngleIndex(angle, angle_combinations_vec);
                    double delta_t = TMath::Abs(time_vec.at(g1) - time_vec.at(g2));

                    // Possible Compton scatter?
                    bool comp_scatter_candidate = ComptonScatterCandidate(angle_index, energy_vec.at(g1), energy_vec.at(g2));

                    // Timing information
                    hist_1D["delta_t"]->Fill(delta_t);

                    // Prompt coincidences
                    if (delta_t < prompt_time_max) {
                        // 1D

                        // 2D
                        hist_2D["sum_energy_angle"]->Fill(angle_index, energy_vec.at(g1) + energy_vec.at(g2));
                        hist_2D["gg_singles"]->Fill(energy_vec.at(g1), energy_vec.at(g2));
                        hist_2D["gg_singles"]->Fill(energy_vec.at(g2), energy_vec.at(g1));
                        hist_2D_prompt[Form("index_%02i_sum", angle_index)]->Fill(energy_vec.at(g1) + energy_vec.at(g2), energy_vec.at(g1));
                        // Compton algorithm
                        if (comp_scatter_candidate){
                            hist_2D_comp_alg_acc[Form("index_%02i", angle_index)]->Fill(energy_vec.at(g1) + energy_vec.at(g2), energy_vec.at(g1));
                        } else {
                            hist_2D_comp_alg_rej[Form("index_%02i", angle_index)]->Fill(energy_vec.at(g1) + energy_vec.at(g2), energy_vec.at(g1));
                        }
                    }
                    if (delta_t > bg_time_min) {
                        // 1D

                        // 2D
                        hist_2D["sum_energy_angle_tr"]->Fill(angle_index, energy_vec.at(g1) + energy_vec.at(g2));
                        hist_2D_tr[Form("index_%02i_sum_tr", angle_index)]->Fill(energy_vec.at(g1) + energy_vec.at(g2), delta_t);
                        if (IsInSlice(delta_t, prompt_time_max)) {
                            // fill histogram and scale by 1/5 since we want the average of 5 time slices
                            hist_2D_tr[Form("index_%02i_sum_tr_avg", angle_index)]->Fill(energy_vec.at(g1) + energy_vec.at(g2), energy_vec.at(g1), 1.0 / 5.0);
                        }
                    }

                } // grif2
            } // grif1
        } // energy_vec


        if (i % 10000 == 0) {
            progress_bar.display();
        }
        ++progress_bar;         // iterates progress_bar


        // cleaning up for next event
        energy_vec.clear();
        pos_vec.clear();
        time_vec.clear();
        kvalue_vec.clear();
        detector_vec.clear();

    }     // end TChain loop
    progress_bar.done();
} // FillHistograms()

/************************************************************//**
 * Pre process data
 *
 ***************************************************************/
void HistogramManager::PreProcessData()
{
    int det_id = -1;
    bool multiplicity_limit_bool = true;
    int multiplicity_limit = 2;

    energy_vec.clear();
    pos_vec.clear();
    time_vec.clear();
    kvalue_vec.clear();
    detector_vec.clear();

    if (multiplicity_limit_bool) {
        if (fGrif->GetSuppressedMultiplicity(fGriffinBgo) == multiplicity_limit) {         // multiplicity filter
            for (auto j = 0; j < fGrif->GetSuppressedMultiplicity(fGriffinBgo); ++j) {
                det_id = fGrif->GetSuppressedHit(j)->GetArrayNumber();
                if (det_id == -1) {
                    std::cout << "BAD DETECTOR" << std::endl;
                    continue;
                }
                //if(fGrif->GetSuppressedHit(j)->GetKValue()!=700) {continue;} // removes GRIFFIN hits pileup events
                float temp_energy = fGrif->GetSuppressedHit(j)->GetEnergy();
                if (temp_energy < 5.0) {
                    continue;
                }

                energy_vec.push_back(temp_energy);
                pos_vec.push_back(fGrif->GetSuppressedHit(j)->GetPosition(145.0));
                time_vec.push_back(fGrif->GetSuppressedHit(j)->GetTime());
                kvalue_vec.push_back(fGrif->GetSuppressedHit(j)->GetKValue());
                detector_vec.push_back(det_id);
            }
        }         // multiplicity filter
    } else {
        for (auto j = 0; j < fGrif->GetSuppressedMultiplicity(fGriffinBgo); ++j) {
            det_id = fGrif->GetSuppressedHit(j)->GetArrayNumber();
            if (det_id == -1) {
                std::cout << "BAD DETECTOR" << std::endl;
                continue;
            }
            //if(fGrif->GetSuppressedHit(j)->GetKValue()!=700) {continue;}     // removes GRIFFIN hits pileup events

            // Applying secondary linear energy calibration
            //energy_temp = offsets[det_id-1] + gains[det_id-1]*fGrif->GetSuppressedHit(j)->GetEnergy();
            // adds small random number to allow proper rebinning
            //energy_temp += ((double) rand() / RAND_MAX - 0.5);

            energy_vec.push_back(fGrif->GetSuppressedHit(j)->GetEnergy());
            pos_vec.push_back(fGrif->GetSuppressedHit(j)->GetPosition(145.0));
            time_vec.push_back(fGrif->GetSuppressedHit(j)->GetTime());
            kvalue_vec.push_back(fGrif->GetSuppressedHit(j)->GetKValue());
            detector_vec.push_back(det_id);
        }
    }
} // PreProcessData

/************************************************************//**
 * Checks if the angle and energy between two photons can
 * be a Compton scatter
 *
 * @param angle Angle between two photons (radians)
 * @param energy_1 Energy of first photon (keV)
 * @param energy_2 Energy of second photon (keV)
 *****************************************************************************/
bool HistogramManager::ComptonScatterCandidate(int angle_index, float energy_1, float energy_2)
{
    float energy_total = energy_1 + energy_2;
    float compton_angle_1 = compton_limit_vec_ftb.at(angle_index);
    float compton_angle_2 = compton_limit_vec_btf.at(angle_index);

    float energy_scatter_high = ComptonScatter(compton_angle_1, energy_total);
    float energy_scatter_low = ComptonScatter(compton_angle_2, energy_total);

    if (((energy_scatter_low < energy_1) && (energy_1 < energy_scatter_high)) || ((energy_scatter_low < energy_2) && (energy_2 < energy_scatter_high))){
        return true;
    } else {
        return false;
    }
} // end ComptonScatterCandidate()

/************************************************************//**
 * Compton scatter formula
 *
 * @param angle Compton scatter angle (rad)
 * @param energy Energy of initial photon (keV)
 *****************************************************************************/
float HistogramManager::ComptonScatter(double angle, float energy)
{
    float electron_rest_mass_energy = 511.; //keV
    return energy / (1 + (energy / electron_rest_mass_energy) * (1 - TMath::Cos(angle)));
} // end ComptonScatter()

/************************************************************//**
 * Returns the angular index
 *
 * @param angle The angle between two gammas
 * @param vec Vector of angles
 *****************************************************************************/
int HistogramManager::GetAngleIndex(double angle, std::vector<double> vec)
{

    // corner cases
    if (angle <= vec.front()) { return 0;}
    if (angle >= vec.back() - 1.) { return vec.size() - 1;}

    // binary search
    unsigned int i = 0, j = vec.size(), mid = 0;
    while ( i < j ) {
        mid = (i + j) / 2;
        if (vec[mid] == angle) return vec[mid];
        // searching left half
        if (angle < vec[mid]) {
            // if angle is greater than previous to mid, return closest of two
            if (mid > 0 && angle > vec[mid - 1]) {
                return GetClosest(mid - 1, mid, angle_combinations_vec, angle);
            }
            // repeat for left half
            j = mid;
        }
        // if angle is greater than mid
        else{
            if (mid < vec.size() - 1 && angle < vec[mid + 1]) {
                return GetClosest(mid, mid + 1, angle_combinations_vec, angle);
            }
            // update i
            i = mid + 1;
        }
    }
    // Only single element left after search
    return mid;
} // GetAngleIndex

/************************************************************//**
 * Returns the value closest to the target
 * Assumes val2 is greater than val1 and target lies inbetween the two
 *
 * @param val1 First value to compare
 * @param val2 Second value to compare
 * @param vec Vector of values
 * @param target Target value
 *****************************************************************************/
int HistogramManager::GetClosest(int val1, int val2, std::vector<double> vec, double target)
{
    if ((target - vec[val1]) >= (vec[val2] - target))
        return val2;
    else
        return val1;
} // GetClosest

/************************************************************//**
 * Opens Root files
 *
 * @param file_name Analysis file path
 ***************************************************************/
void HistogramManager::WriteHistogramsToFile()
{
    TFile *out_file = new TFile("histograms.root", "RECREATE");
    TDirectory *time_random_dir = out_file->mkdir("time_random");
    TDirectory *prompt_angle_dir = out_file->mkdir("prompt_angle");
    TDirectory *comp_dir_rej = out_file->mkdir("compton_alg_rej");
    TDirectory *comp_dir_acc = out_file->mkdir("compton_alg_acc");
    std::cout << "Writing histograms to file: " << out_file->GetName() << std::endl;

    out_file->cd();
    for(auto my_histogram : hist_1D) {
        my_histogram.second->Write();
    }
    for(auto my_histogram : hist_2D) {
        my_histogram.second->Write();
    }
    time_random_dir->cd();
    for(auto my_histogram : hist_2D_tr) {
        my_histogram.second->Write();
    }
    prompt_angle_dir->cd();
    for(auto my_histogram : hist_2D_prompt) {
        my_histogram.second->Write();
    }
    comp_dir_acc->cd();
    for(auto my_histogram : hist_2D_comp_alg_acc) {
        my_histogram.second->Write();
    }
    comp_dir_rej->cd();
    for(auto my_histogram : hist_2D_comp_alg_rej) {
        my_histogram.second->Write();
    }
    out_file->Close();
} // WriteHistogramsToFile()

/************************************************************//**
 * Checks if time difference in is background slices
 *
 ***************************************************************/
bool HistogramManager::IsInSlice(double delta_t, int prompt_time)
{
    int slice_edges[5] = {510, 617, 725, 832, 940};

    for (auto i = 0; i < 5; i++) {
        if ((delta_t > slice_edges[i]) && (delta_t < slice_edges[i] + prompt_time)) {
            return true;
        } else {
            continue;
        }
    }

    // if the time difference does not fall in the bg slices
    return false;

} // end IsInSlice()
