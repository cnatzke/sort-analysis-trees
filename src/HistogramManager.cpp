//////////////////////////////////////////////////////////////////////////////////
// Creates and fills histograms
//
// Author:          Connor Natzke (cnatzke@triumf.ca)
// Creation Date: Wednesday July 29, 2020	T15:22:33-07:00
// Last Update:   Wednesday July 29, 2020	T15:22:45-07:00
// Usage:
//
//////////////////////////////////////////////////////////////////////////////////
#include <map>
#include <algorithm>
#include "HistogramManager.h"
#include "progress_bar.h"
#include "LoadingMessenger.h"

/************************************************************//**
 * Constructor
 ***************************************************************/
HistogramManager::HistogramManager(std::string compton_limits_file) : _compton_limits_file(compton_limits_file)
{
} // end Constructor

/************************************************************//**
 * Destructor
 ***************************************************************/
HistogramManager::~HistogramManager(void)
{
} // end Destructor

/************************************************************//**
 * Creates and Fills histograms
 *
 * @param verbose Verbosity level
 ***************************************************************/
void HistogramManager::MakeHistograms(TChain *input_chain)
{
    int verbosity = 1;
    _detector_radius = 145.0; //mm
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

    hist_1D["singles_energy"] = new TH1D("singles_energy", "gamma singles", energy_bins_max, energy_bins_min, energy_bins_max);
    hist_1D["addback_energy"] = new TH1D("addback_energy", "Addback Energy", energy_bins_max, energy_bins_min, energy_bins_max);
    hist_1D["singles_unsup_energy"] = new TH1D("singles_unsup_energy", "Unsuppressed singles energy", energy_bins_max, energy_bins_min, energy_bins_max);
    hist_1D["singles_reconstructed_energy"] = new TH1D("singles_reconstructed_energy", "singles Energy", energy_bins_max, energy_bins_min, energy_bins_max);
    hist_1D["addback_reconstructed_energy"] = new TH1D("addback_reconstructed_energy", "Addback Energy", energy_bins_max, energy_bins_min, energy_bins_max);
    hist_1D["singles_rejected_energy"] = new TH1D("singles_rejected_energy", "singles Energy", energy_bins_max, energy_bins_min, energy_bins_max);
    hist_1D["addback_rejected_energy"] = new TH1D("addback_rejected_energy", "Addback Energy", energy_bins_max, energy_bins_min, energy_bins_max);
    // accepted Compton energies
    hist_1D["singles_compton_energy"] = new TH1D("singles_compton_energy", "Accepted Compton gamma singles", energy_bins_max, energy_bins_min, energy_bins_max);
    hist_1D["addback_compton_energy"] = new TH1D("addback_compton_energy", "Accepted Compton addback energy", energy_bins_max, energy_bins_min, energy_bins_max);


    // 2D Histograms
    if (verbose > 0 ) std::cout << "Creating 2D histograms ... " << std::endl;
    hist_2D["singles_energy_channel"] = new TH2D("singles_energy_channel", "#gamma singles vs channel", 70, 0, 70, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["singles_gg_matrix"] = new TH2D("singles_gg_matrix", "", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["singles_reconstructed_gg_matrix"] = new TH2D("singles_reconstructed_gg_matrix", "", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["singles_rejected_gg_matrix"] = new TH2D("singles_rejected_gg_matrix", "", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["singles_unsup_gg_matrix"] = new TH2D("singles_unsup_gg_matrix", "", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["singles_sum_angle"] = new TH2D("singles_sum_angle", "Sum Energy vs angle index", 70, 0, 70, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["sum_energy_angle_tr"] = new TH2D("sum_energy_angle_tr", "Sum Energy vs angle index (time random);angle index;sum energy [keV]", 70, 0, 70, energy_bins_max, energy_bins_min, energy_bins_max);

    hist_2D["addback_energy_channel"] = new TH2D("addback_energy_channel", "#gamma addback vs channel", 70, 0, 70, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["addback_gg_matrix"] = new TH2D("addback_gg_matrix", "", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["addback_sum_angle"] = new TH2D("addback_sum_angle", "Sum Energy vs angle index", 70, 0, 70, energy_bins_max, energy_bins_min, energy_bins_max);

    hist_2D["addback_reconstructed_energy_channel"] = new TH2D("addback_reconstructed_energy_channel", "#gamma addback_reconstructed vs channel", 70, 0, 70, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["addback_reconstructed_gg_matrix"] = new TH2D("addback_reconstructed_gg_matrix", "", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["addback_reconstructed_sum_angle"] = new TH2D("addback_reconstructed_sum_angle", "Sum Energy vs angle index", 70, 0, 70, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["addback_rejected_gg_matrix"] = new TH2D("addback_rejected_gg_matrix", "", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);

    hist_2D["singles_compton_gg_matrix"] = new TH2D("singles_compton_gg_matrix", "", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);
    hist_2D["addback_compton_gg_matrix"] = new TH2D("addback_compton_gg_matrix", "", energy_bins_max, energy_bins_min, energy_bins_max, energy_bins_max, energy_bins_min, energy_bins_max);


} // InitializeHistograms()

/************************************************************//**
 * Fills histograms
 *
 * @param gChain Data chain
 ***************************************************************/
void HistogramManager::FillHistograms(TChain *gChain)
{

    _prompt_time_max = 30; // ns
    _bg_time_min = 500; // ns

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

    ComptonRecovery * comp_check = new ComptonRecovery(_compton_limits_file);

    // display funny loading message
    LoadingMessenger load_man;
    load_man.DisplayLoadingMessage();

    long analysis_entries = gChain->GetEntries();
    ProgressBar progress_bar(analysis_entries, 70, '=', ' ');
    //for (auto i = 0; i < 100; i++) {
    for (auto i = 0; i < analysis_entries; i++) {
        gChain->GetEntry(i);
        _event_number = i;

        // Applies multiplicity filter and recovers intra-clover compton scatters
        PreProcessData(comp_check);

        // suppressed singles
        if (singles_energy_vec.size() > 0) {
            for (auto g1 = 0; g1 < (int) singles_energy_vec.size(); g1++) {
                hist_1D["k_value"]->Fill(singles_kvalue_vec.at(g1));
                hist_1D["singles_energy"]->Fill(singles_energy_vec.at(g1));
                hist_2D["singles_energy_channel"]->Fill(singles_id_vec.at(g1), singles_energy_vec.at(g1));

                //for (auto g2 = 0; g2 < singles_energy_vec.size(); g2++) { // symmetric matrices
                for (auto g2 = g1 + 1; g2 < (int) singles_energy_vec.size(); g2++) {     // asymmetric looping
                    if (g1 == g2) continue;

                    double angle = singles_pos_vec.at(g1).Angle(singles_pos_vec.at(g2));
                    int angle_index = GetAngleIndex(angle * rad_to_degree, angle_combinations_vec);
                    if (angle * rad_to_degree < 0.0001 || angle * rad_to_degree > 180.0) {
                        continue;
                    }

                    double delta_t = TMath::Abs(singles_time_vec.at(g1) - singles_time_vec.at(g2));
                    // Prompt coincidences
                    if (delta_t < _prompt_time_max) {
                        // 1D

                        // 2D
                        hist_2D["singles_sum_angle"]->Fill(angle_index, singles_energy_vec.at(g1) + singles_energy_vec.at(g2));
                        hist_2D["singles_gg_matrix"]->Fill(singles_energy_vec.at(g1), singles_energy_vec.at(g2));
                        hist_2D["singles_gg_matrix"]->Fill(singles_energy_vec.at(g2), singles_energy_vec.at(g1));
                    } // end prompt coincidence
                } // end g2
            } // end g1
        } // end singles

        // addback
        if (addback_energy_vec.size() > 0) {
            for (auto g1 = 0; g1 < (int) addback_energy_vec.size(); g1++) {
                hist_1D["addback_energy"]->Fill(addback_energy_vec.at(g1));
                hist_2D["addback_energy_channel"]->Fill(addback_id_vec.at(g1), addback_energy_vec.at(g1));

                //for (auto g2 = 0; g2 < addback_energy_vec.size(); g2++) { // symmetric matrices
                for (auto g2 = g1 + 1; g2 < (int) addback_energy_vec.size(); g2++) {     // asymmetric looping
                    if (g1 == g2) continue;

                    double angle = addback_pos_vec.at(g1).Angle(addback_pos_vec.at(g2));
                    int angle_index = GetAngleIndex(angle * rad_to_degree, angle_combinations_vec);
                    if (angle * rad_to_degree < 0.0001 || angle * rad_to_degree > 180.0) {
                        continue;
                    }

                    double delta_t = TMath::Abs(addback_time_vec.at(g1) - addback_time_vec.at(g2));
                    // Prompt coincidences
                    if (delta_t < _prompt_time_max) {
                        // 1D

                        // 2D
                        hist_2D["addback_sum_angle"]->Fill(angle_index, addback_energy_vec.at(g1) + addback_energy_vec.at(g2));
                        hist_2D["addback_gg_matrix"]->Fill(addback_energy_vec.at(g1), addback_energy_vec.at(g2));
                        hist_2D["addback_gg_matrix"]->Fill(addback_energy_vec.at(g2), addback_energy_vec.at(g1));
                    } // end prompt coincidence
                } // end g2
            } // end g1
        } // end addback

        // Compton recovered singles
        if (singles_reconstructed_energy_vec.size() > 0) {
            for (auto g1 = 0; g1 < (int) singles_reconstructed_energy_vec.size(); g1++) {
                hist_1D["singles_reconstructed_energy"]->Fill(singles_reconstructed_energy_vec.at(g1));
                for (auto g2 = g1 + 1; g2 < (int) singles_reconstructed_energy_vec.size(); g2++) {     // asymmetric looping
                    if (g1 == g2) continue;

                    double angle = singles_reconstructed_pos_vec.at(g1).Angle(singles_reconstructed_pos_vec.at(g2));
                    int angle_index = GetAngleIndex(angle * rad_to_degree, angle_combinations_vec);
                    if (angle * rad_to_degree < 0.0001 || angle * rad_to_degree > 180.0) {
                        continue;
                    }

                    double delta_t = TMath::Abs(singles_reconstructed_time_vec.at(g1) - singles_reconstructed_time_vec.at(g2));
                    // Prompt coincidences
                    if (delta_t < _prompt_time_max) {
                        // 1D

                        // 2D
                        hist_2D["singles_reconstructed_gg_matrix"]->Fill(singles_reconstructed_energy_vec.at(g1), singles_reconstructed_energy_vec.at(g2));
                        hist_2D["singles_reconstructed_gg_matrix"]->Fill(singles_reconstructed_energy_vec.at(g2), singles_reconstructed_energy_vec.at(g1));
                    } // end prompt coincidence
                } // end g2
            } // end g1
        } // end compton reconstructed singles

        // Compton recovered addback
        if (addback_reconstructed_energy_vec.size() > 0) {
            for (auto g1 = 0; g1 < (int) addback_reconstructed_energy_vec.size(); g1++) {
                hist_1D["addback_reconstructed_energy"]->Fill(addback_reconstructed_energy_vec.at(g1));
                hist_2D["addback_reconstructed_energy_channel"]->Fill(addback_reconstructed_id_vec.at(g1), addback_reconstructed_energy_vec.at(g1));

                //for (auto g2 = 0; g2 < addback_reconstructed_energy_vec.size(); g2++) { // symmetric matrices
                for (auto g2 = g1 + 1; g2 < (int) addback_reconstructed_energy_vec.size(); g2++) {     // asymmetric looping
                    if (g1 == g2) continue;

                    double angle = addback_reconstructed_pos_vec.at(g1).Angle(addback_reconstructed_pos_vec.at(g2));
                    int angle_index = GetAngleIndex(angle * rad_to_degree, angle_combinations_vec);
                    if (angle * rad_to_degree < 0.0001 || angle * rad_to_degree > 180.0) {
                        continue;
                    }

                    double delta_t = TMath::Abs(addback_reconstructed_time_vec.at(g1) - addback_reconstructed_time_vec.at(g2));
                    // Prompt coincidences
                    if (delta_t < _prompt_time_max) {
                        // 1D

                        // 2D
                        hist_2D["addback_reconstructed_sum_angle"]->Fill(angle_index, addback_reconstructed_energy_vec.at(g1) + addback_reconstructed_energy_vec.at(g2));
                        hist_2D["addback_reconstructed_gg_matrix"]->Fill(addback_reconstructed_energy_vec.at(g1), addback_reconstructed_energy_vec.at(g2));
                        hist_2D["addback_reconstructed_gg_matrix"]->Fill(addback_reconstructed_energy_vec.at(g2), addback_reconstructed_energy_vec.at(g1));
                    } // end prompt coincidence
                } // end g2
            } // end g1
        } // end compton recovered addback

        // Unsuppressed singles
        if (singles_unsup_energy_vec.size() > 0) {
            for (auto g1 = 0; g1 < (int) singles_unsup_energy_vec.size(); g1++) {
                hist_1D["singles_unsup_energy"]->Fill(singles_unsup_energy_vec.at(g1));
                for (auto g2 = g1 + 1; g2 < (int) singles_unsup_energy_vec.size(); g2++) { // asymmetric looping
                    if (g1 == g2) continue;

                    double delta_t = TMath::Abs(singles_unsup_time_vec.at(g1) - singles_unsup_time_vec.at(g2));
                    // Prompt coincidences
                    if (delta_t < _prompt_time_max) {
                        // 1D

                        // 2D
                        hist_2D["singles_unsup_gg_matrix"]->Fill(singles_unsup_energy_vec.at(g1), singles_unsup_energy_vec.at(g2));
                        hist_2D["singles_unsup_gg_matrix"]->Fill(singles_unsup_energy_vec.at(g2), singles_unsup_energy_vec.at(g1));
                    } // end prompt coincidence
                } // end g2
            } // end g1
        } // end compton reconstructed singles

        // singles Compton events
        if (singles_compton_energy_vec.size() > 0) {
            for (auto g1 = 0; g1 < (int) singles_compton_energy_vec.size(); g1++) {
                hist_1D["singles_compton_energy"]->Fill(singles_compton_energy_vec.at(g1));
                for (auto g2 = g1 + 1; g2 < (int) singles_compton_energy_vec.size(); g2++) { // asymmetric looping
                    if (g1 == g2) continue;

                    double delta_t = TMath::Abs(singles_compton_time_vec.at(g1) - singles_compton_time_vec.at(g2));
                    // Prompt coincidences
                    if (delta_t < _prompt_time_max) {
                        // 1D
                        // 2D
                        hist_2D["singles_compton_gg_matrix"]->Fill(singles_compton_energy_vec.at(g1), singles_compton_energy_vec.at(g2));
                        hist_2D["singles_compton_gg_matrix"]->Fill(singles_compton_energy_vec.at(g2), singles_compton_energy_vec.at(g1));
                    } // end prompt coincidence
                } // end g2
            } // end g1
        } // end singles Compton events

        // addback Compton events
        if (addback_compton_energy_vec.size() > 0) {
            for (auto g1 = 0; g1 < (int) addback_compton_energy_vec.size(); g1++) {
                hist_1D["addback_compton_energy"]->Fill(addback_compton_energy_vec.at(g1));
                for (auto g2 = g1 + 1; g2 < (int) addback_compton_energy_vec.size(); g2++) { // asymmetric looping
                    if (g1 == g2) continue;

                    double delta_t = TMath::Abs(addback_compton_time_vec.at(g1) - addback_compton_time_vec.at(g2));
                    // Prompt coincidences
                    if (delta_t < _prompt_time_max) {
                        // 1D
                        // 2D
                        hist_2D["addback_compton_gg_matrix"]->Fill(addback_compton_energy_vec.at(g1), addback_compton_energy_vec.at(g2));
                        hist_2D["addback_compton_gg_matrix"]->Fill(addback_compton_energy_vec.at(g2), addback_compton_energy_vec.at(g1));
                    } // end prompt coincidence
                } // end g2
            } // end g1
        } // end addback Compton events

        // singles Compton rejected events
        if (singles_rejected_energy_vec.size() > 0) {
            for (auto g1 = 0; g1 < (int) singles_rejected_energy_vec.size(); g1++) {
                hist_1D["singles_rejected_energy"]->Fill(singles_rejected_energy_vec.at(g1));
                for (auto g2 = g1 + 1; g2 < (int) singles_rejected_energy_vec.size(); g2++) { // asymmetric looping
                    if (g1 == g2) continue;

                    double delta_t = TMath::Abs(singles_rejected_time_vec.at(g1) - singles_rejected_time_vec.at(g2));
                    // Prompt coincidences
                    if (delta_t < _prompt_time_max) {
                        // 1D
                        // 2D
                        hist_2D["singles_rejected_gg_matrix"]->Fill(singles_rejected_energy_vec.at(g1), singles_rejected_energy_vec.at(g2));
                        hist_2D["singles_rejected_gg_matrix"]->Fill(singles_rejected_energy_vec.at(g2), singles_rejected_energy_vec.at(g1));
                    } // end prompt coincidence
                } // end g2
            } // end g1
        } // end singles rejected events

        // addback Compton rejected events
        if (addback_rejected_energy_vec.size() > 0) {
            for (auto g1 = 0; g1 < (int) addback_rejected_energy_vec.size(); g1++) {
                hist_1D["addback_rejected_energy"]->Fill(addback_rejected_energy_vec.at(g1));
                for (auto g2 = g1 + 1; g2 < (int) addback_rejected_energy_vec.size(); g2++) { // asymmetric looping
                    if (g1 == g2) continue;

                    double delta_t = TMath::Abs(addback_rejected_time_vec.at(g1) - addback_rejected_time_vec.at(g2));
                    // Prompt coincidences
                    if (delta_t < _prompt_time_max) {
                        // 1D
                        // 2D
                        hist_2D["addback_rejected_gg_matrix"]->Fill(addback_rejected_energy_vec.at(g1), addback_rejected_energy_vec.at(g2));
                        hist_2D["addback_rejected_gg_matrix"]->Fill(addback_rejected_energy_vec.at(g2), addback_rejected_energy_vec.at(g1));
                    } // end prompt coincidence
                } // end g2
            } // end g1
        } // end addback rejected events

        if (i % 10000 == 0) {
            progress_bar.display();
        }
        ++progress_bar;         // iterates progress_bar


    }     // end TChain loop
    progress_bar.done();

    delete comp_check;
} // FillHistograms()



/************************************************************//**
 * Pre process data
 *
 ***************************************************************/
void HistogramManager::PreProcessData(ComptonRecovery * comp_check)
{
    int det_id = -1;
    bool diagnostic_verbosity = false;

    // energy vectors
    singles_energy_vec.clear();
    addback_energy_vec.clear();
    singles_unsup_energy_vec.clear();
    singles_reconstructed_energy_vec.clear();
    addback_reconstructed_energy_vec.clear();
    singles_rejected_energy_vec.clear();
    addback_rejected_energy_vec.clear();
    singles_compton_energy_vec.clear();
    addback_compton_energy_vec.clear();
    // position vectors
    singles_pos_vec.clear();
    addback_pos_vec.clear();
    singles_unsup_pos_vec.clear();
    singles_reconstructed_pos_vec.clear();
    addback_reconstructed_pos_vec.clear();
    // time vectors
    singles_time_vec.clear();
    addback_time_vec.clear();
    singles_unsup_time_vec.clear();
    singles_reconstructed_time_vec.clear();
    addback_reconstructed_time_vec.clear();
    singles_rejected_time_vec.clear();
    addback_rejected_time_vec.clear();
    singles_compton_time_vec.clear();
    addback_compton_time_vec.clear();
    // detector id vectors
    singles_id_vec.clear();
    addback_id_vec.clear();
    singles_clover_id_vec.clear();
    singles_unsup_id_vec.clear();
    singles_reconstructed_id_vec.clear();
    addback_reconstructed_id_vec.clear();
    // misc vectors
    singles_kvalue_vec.clear();

    // unsuppressed singles
    for(auto g1 = 0; g1 < fGrif->GetMultiplicity(); g1++) {
        TGriffinHit * grif1 = static_cast<TGriffinHit*>(fGrif->GetHit(g1));

        singles_unsup_energy_vec.push_back(grif1->GetEnergy());
        singles_unsup_pos_vec.push_back(grif1->GetPosition(_detector_radius));
        singles_unsup_time_vec.push_back(grif1->GetTime());
        singles_unsup_id_vec.push_back(grif1->GetArrayNumber());
    } // end unsuppressed singles

    // suppressed singles
    int true_singles_multiplicity = static_cast<int>(fGrif->GetSuppressedMultiplicity(fGriffinBgo));
    std::vector<int> accepted_singles_compton_indices;
    bool found_singles_reconstruction_event = false;
    for(auto g1 = 0; g1 < true_singles_multiplicity; g1++) {
        TGriffinHit * grif1 = static_cast<TGriffinHit*>(fGrif->GetSuppressedHit(g1));

        singles_energy_vec.push_back(grif1->GetEnergy());
        singles_pos_vec.push_back(grif1->GetPosition(_detector_radius));
        singles_time_vec.push_back(grif1->GetTime());
        singles_id_vec.push_back(grif1->GetArrayNumber());
        singles_clover_id_vec.push_back(grif1->GetDetector());
        singles_kvalue_vec.push_back(grif1->GetKValue());

        // Compton recovery logic
        for(auto g2 = g1 + 1; g2 < true_singles_multiplicity; g2++) { // loop MUST be assymmetric because of Compton recovery alg.

            // only check for reconstructed indices if one has been found
            if (diagnostic_verbosity && found_singles_reconstruction_event && true_singles_multiplicity > 3) std::cout << _event_number << " | " << g1 << " | " << g2 << " | "<< std::endl;
            if (found_singles_reconstruction_event) {
                if (std::find(accepted_singles_compton_indices.begin(), accepted_singles_compton_indices.end(), g1) != accepted_singles_compton_indices.end()) {
                    continue;
                }
                if (std::find(accepted_singles_compton_indices.begin(), accepted_singles_compton_indices.end(), g2) != accepted_singles_compton_indices.end()) {
                    continue;
                }
            }
            if (diagnostic_verbosity && found_singles_reconstruction_event && true_singles_multiplicity > 3) std::cout << _event_number << " | " << g1 << " | " << g2 << " | Passed"<< std::endl;

            TGriffinHit * grif2 = static_cast<TGriffinHit*>(fGrif->GetSuppressedHit(g2));

            // find angle between hits
            double angle = grif1->GetPosition(_detector_radius).Angle(grif2->GetPosition(_detector_radius));
            int angle_index = GetAngleIndex(angle * rad_to_degree, angle_combinations_vec);
            // angle safety check
            if (angle * rad_to_degree < 0.0001 || angle * rad_to_degree > 180.0) {
                continue;
            }

            double delta_t = TMath::Abs(grif1->GetTime() - grif2->GetTime());
            if (delta_t < _prompt_time_max) {
                // check for possible intra-clover Compton scatter
                bool compton_scatter_candidate = comp_check->ComptonScatterCandidate(angle_index, grif1->GetEnergy(), grif2->GetEnergy());
                if (compton_scatter_candidate) {
                    singles_reconstructed_energy_vec.push_back(grif1->GetEnergy() + grif2->GetEnergy());
                    singles_compton_energy_vec.push_back(grif1->GetEnergy());
                    singles_compton_time_vec.push_back(grif1->GetTime());
                    singles_compton_energy_vec.push_back(grif2->GetEnergy());
                    singles_compton_time_vec.push_back(grif2->GetTime());
                    // assign reconstructed energy to clover with highest energy
                    if (comp_check->FirstHitHigh(grif1, grif2)) {
                        singles_reconstructed_pos_vec.push_back(grif1->GetPosition(_detector_radius));
                        singles_reconstructed_id_vec.push_back(grif1->GetArrayNumber());
                    } else{
                        singles_reconstructed_pos_vec.push_back(grif2->GetPosition(_detector_radius));
                        singles_reconstructed_id_vec.push_back(grif2->GetArrayNumber());
                    }
                    // assign time as time of first hit
                    singles_reconstructed_time_vec.push_back(comp_check->GetReconstructedTime(grif1, grif2));

                    // removing hits from future loops
                    accepted_singles_compton_indices.push_back(g1);
                    accepted_singles_compton_indices.push_back(g2);
                    found_singles_reconstruction_event = true;

                    // checking correct events are passed
                    if (diagnostic_verbosity && true_singles_multiplicity > 3) {
                        std::cout << "---> " << _event_number
                                  << " | " << g1
                                  << " | " << g2
                                  << " | " << true_singles_multiplicity
                            //<< " | " << effective_singles_multiplicity
                                  << " | " << std::endl;
                    }
                } // end compton_scatter_candidate
            } // end prompt coincidence
        } // end grif2 & compton recovery logic

        // add first grif hit if it was not a compton candidate with any other hits in the event
        if (std::find(accepted_singles_compton_indices.begin(), accepted_singles_compton_indices.end(), g1) != accepted_singles_compton_indices.end()) {
            continue;
        } else {
            singles_reconstructed_energy_vec.push_back(grif1->GetEnergy());
            singles_reconstructed_pos_vec.push_back(grif1->GetPosition(_detector_radius));
            singles_reconstructed_time_vec.push_back(grif1->GetTime());
            singles_reconstructed_id_vec.push_back(grif1->GetArrayNumber());

            singles_rejected_energy_vec.push_back(grif1->GetEnergy());
            singles_rejected_time_vec.push_back(grif1->GetTime());

            if (diagnostic_verbosity && true_singles_multiplicity > 3) {
                std::cout << _event_number << "| Added hit " << g1 << " to list" << std::endl;
            }
        } // end grif1 reconstruction check
    } // end g1 & singles


    // suppressed addback
    int true_addback_multiplicity = static_cast<int>(fGrif->GetSuppressedAddbackMultiplicity(fGriffinBgo));
    std::vector<int> accepted_addback_compton_indices;
    bool found_reconstruction_event = false;
    for (auto g1 = 0; g1 < true_addback_multiplicity; g1++) {
        TGriffinHit * grif1 = static_cast<TGriffinHit*>(fGrif->GetSuppressedAddbackHit(g1));

        addback_energy_vec.push_back(grif1->GetEnergy());
        addback_pos_vec.push_back(grif1->GetPosition(_detector_radius));
        addback_time_vec.push_back(grif1->GetTime());
        addback_id_vec.push_back(grif1->GetArrayNumber());

        // Compton recovery logic
        for(auto g2 = g1 + 1; g2 < true_addback_multiplicity; g2++) { // loop MUST be assymmetric because of Compton recovery alg.

            // only check for reconstructed indices if one has been found
            if (diagnostic_verbosity && found_reconstruction_event && true_addback_multiplicity > 3) std::cout << _event_number << " | " << g1 << " | " << g2 << " | "<< std::endl;
            if (found_reconstruction_event) {
                if (std::find(accepted_addback_compton_indices.begin(), accepted_addback_compton_indices.end(), g1) != accepted_addback_compton_indices.end()) {
                    continue;
                }
                if (std::find(accepted_addback_compton_indices.begin(), accepted_addback_compton_indices.end(), g2) != accepted_addback_compton_indices.end()) {
                    continue;
                }
            }
            if (diagnostic_verbosity && found_reconstruction_event && true_addback_multiplicity > 3) std::cout << _event_number << " | " << g1 << " | " << g2 << " | Passed"<< std::endl;

            TGriffinHit * grif2 = static_cast<TGriffinHit*>(fGrif->GetSuppressedAddbackHit(g2));

            // find angle between hits
            double angle = grif1->GetPosition(_detector_radius).Angle(grif2->GetPosition(_detector_radius));
            int angle_index = GetAngleIndex(angle * rad_to_degree, angle_combinations_vec);
            // angle safety check
            if (angle * rad_to_degree < 0.0001 || angle * rad_to_degree > 180.0) {
                continue;
            }

            double delta_t = TMath::Abs(grif1->GetTime() - grif2->GetTime());
            if (delta_t < _prompt_time_max) {
                // check for possible intra-clover Compton scatter
                bool compton_scatter_candidate = comp_check->ComptonScatterCandidate(angle_index, grif1->GetEnergy(), grif2->GetEnergy());
                if (compton_scatter_candidate) {
                    addback_reconstructed_energy_vec.push_back(grif1->GetEnergy() + grif2->GetEnergy());
                    addback_compton_energy_vec.push_back(grif1->GetEnergy());
                    addback_compton_time_vec.push_back(grif1->GetTime());
                    addback_compton_energy_vec.push_back(grif2->GetEnergy());
                    addback_compton_time_vec.push_back(grif2->GetTime());
                    // assign reconstructed energy to clover with highest energy
                    if (comp_check->FirstHitHigh(grif1, grif2)) {
                        addback_reconstructed_pos_vec.push_back(grif1->GetPosition(_detector_radius));
                        addback_reconstructed_id_vec.push_back(grif1->GetArrayNumber());
                    } else{
                        addback_reconstructed_pos_vec.push_back(grif2->GetPosition(_detector_radius));
                        addback_reconstructed_id_vec.push_back(grif2->GetArrayNumber());
                    }
                    // assign time as time of first hit
                    addback_reconstructed_time_vec.push_back(comp_check->GetReconstructedTime(grif1, grif2));

                    // removing hits from future loops
                    accepted_addback_compton_indices.push_back(g1);
                    accepted_addback_compton_indices.push_back(g2);
                    found_reconstruction_event = true;

                    // checking correct events are passed
                    if (diagnostic_verbosity && true_addback_multiplicity > 3) {
                        std::cout << "---> " << _event_number
                                  << " | " << g1
                                  << " | " << g2
                                  << " | " << true_addback_multiplicity
                            //<< " | " << effective_addback_multiplicity
                                  << " | " << std::endl;
                    }
                } // end compton_scatter_candidate
            } // end prompt coincidence
        } // end grif2 & compton recovery logic

        // add first grif hit if it was not a compton candidate with any other hits in the event
        if (std::find(accepted_addback_compton_indices.begin(), accepted_addback_compton_indices.end(), g1) != accepted_addback_compton_indices.end()) {
            continue;
        } else {
            addback_reconstructed_energy_vec.push_back(grif1->GetEnergy());
            addback_reconstructed_pos_vec.push_back(grif1->GetPosition(_detector_radius));
            addback_reconstructed_time_vec.push_back(grif1->GetTime());
            addback_reconstructed_id_vec.push_back(grif1->GetArrayNumber());

            addback_rejected_energy_vec.push_back(grif1->GetEnergy());
            addback_rejected_time_vec.push_back(grif1->GetTime());

            if (diagnostic_verbosity && true_addback_multiplicity > 3) {
                std::cout << _event_number << "| Added hit " << g1 << " to list" << std::endl;
            }
        } // end grif1 reconstruction check
    } // end grif1 - addback

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
    std::cout << "Writing histograms to file: " << out_file->GetName() << std::endl;

    out_file->cd();
    for(auto my_histogram : hist_1D) {
        my_histogram.second->Write();
    }
    for(auto my_histogram : hist_2D) {
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
