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
    // create histograms
    InitializeHistograms(verbosity);
    FillHistograms(input_chain);

} // GenerateHistogramFile()

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
    hist_2D["sum_energy_angle_tr"] = new TH2D("sum_energy_angle_tr", "Sum Energy vs angle index (time random)", 70, 0, 70, energy_bins_max, energy_bins_min, energy_bins_max);

    // individual histograms for each angular bin
    for (unsigned int i = 0; i < angle_combinations_vec.size(); i++) {
        hist_2D_prompt[Form("index_%02i_sum", i)] = new TH2D(Form("index_%02i_sum", i), ";sum energy [keV];#gamma_1 energy [keV]", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
        hist_2D_tr[Form("index_%02i_sum_tr", i)] = new TH2D(Form("index_%02i_sum_tr", i), ";sum energy [keV];#Deltat [ns]", energy_bins_max, energy_bins_min, energy_bins_max, 600, 400, 1000);
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

                    double angle = pos_vec.at(g1).Angle(pos_vec.at(g2)) * rad_to_degree;
                    if (angle < 0.0001 || angle > 180.0) {
                        continue;
                    }

                    int angle_index = GetAngleIndex(angle, angle_combinations_vec);
                    double delta_t = TMath::Abs(time_vec.at(g1) - time_vec.at(g2));

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
                    }
                    if (delta_t > bg_time_min) {
                        // 1D

                        // 2D
                        hist_2D["sum_energy_angle_tr"]->Fill(angle_index, energy_vec.at(g1) + energy_vec.at(g2));
                        hist_2D_tr[Form("index_%02i_sum_tr", angle_index)]->Fill(energy_vec.at(g1) + energy_vec.at(g2), delta_t);
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
 * Applies linear calibration to data points
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
    TDirectory *time_random_dir = out_file->mkdir("time-random");
    TDirectory *prompt_angle_dir = out_file->mkdir("prompt-angle");
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
    out_file->Close();
} // WriteHistogramsToFile()
