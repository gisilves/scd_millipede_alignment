#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <limits>
#include <cmath>
#include "TLinearFitter.h"
#include "Mille.h"
#include "TFile.h"
#include "TH1.h"

const std::vector<double> Z_LAYERS = {0.0, 5.0, 55.0, 60.0, 110.0, 115.0, 165.0, 170.0};
const std::vector<bool> IS_X_LAYER = {false, true, false, true, false, true, false, true};
const std::vector<bool> IS_ACTIVE = {true, false, true, true, true, true, false, true};
const float SIGMA_HIT = 0.05;

double cY_ref = 0.085;   // main-peak position from Residual_0
double cX_ref = -0.04;    // read this off Residual_3 (or 5/7) for the X view
const double CUT    = 0.05;   // several core widths; widen or tighten as needed

const int GLOBAL_LABELS[8][4] ={
    {1012, 1013, 1014, 1015}, // Layer 0 (Y, inactive)
    {2011, 2013, 2014, 2015}, // Layer 1 (X, active)
    {3012, 3013, 3014, 3015}, // Layer 2 (Y, active)
    {4011, 4013, 4014, 4015}, // Layer 3 (X, active)
    {5012, 5013, 5014, 5015}, // Layer 4 (Y, active)
    {6011, 6013, 6014, 6015}, // Layer 5 (X, active)
    {7012, 7013, 7014, 7015}, // Layer 6 (Y, inactive)
    {8011, 8013, 8014, 8015}  // Layer 7 (X, active)
};

std::map<int, double> alignment_offsets;

void loadGeometry()
{
    std::ifstream file("current_geometry.txt");
    if (!file.is_open())
    {
        std::cout << "--> Missing current_geometry.txt file. Using nominal positions (Iteration 1)." << std::endl;
        return;
    }
    int id;
    double offset;
    while (file >> id >> offset)
    {
        alignment_offsets[id] = offset;
    }
    file.close();
    std::cout << "--> Geometry loaded successfully." << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <iteration> hits_file" << std::endl;
        return 1;
    }

    int iteration = (argc > 1) ? std::atoi(argv[1]) : 1;
    std::string hits_file = (argc > 2) ? argv[2] : "hits_misaligned.csv";

    // 1. Load the alignment corrections calculated in previous steps
    loadGeometry();
    std::cout << "Loaded geometry" << std::endl;

    std::ifstream csvFile(hits_file);
    if (!csvFile.is_open())
    {
        std::cerr << "Error: cannot open the CSV file with the measurements." << std::endl;
        return 1;
    }

    // 2. Output file for residual histograms

    // Rootfiles are named as "residuals_<iteration>.root"
    std::string output_file = "rootfiles/residuals_" + std::to_string(iteration) + ".root";
    TFile outputFile(output_file.c_str(), "RECREATE");
    std::vector<TH1*> residuals;
    for (int i = 0; i < 8; ++i)
    {
        std::string name = "Residual_" + std::to_string(i);
        residuals.push_back(new TH1F(name.c_str(), name.c_str(), 2500, -1.0, 1.0));
    }

    Mille milleWriter("mille_output_run001.bin");
    std::string line;
    int eventCount = 0;

    TLinearFitter fitterX(1, "pol1");
    TLinearFitter fitterY(1, "pol1");

    std::vector<int> idxX = {1, 3, 5, 7};
    std::vector<int> idxY = {0, 2, 4, 6};

    while (std::getline(csvFile, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::string val;
        std::vector<double> raw_hits;

        while (std::getline(ss, val, ','))
        {
            // Trim leading/trailing whitespace
            size_t start = val.find_first_not_of(" \t\r\n");
            size_t end = val.find_last_not_of(" \t\r\n");
            
            if (start == std::string::npos) {
                // Empty field
                raw_hits.push_back(std::numeric_limits<double>::quiet_NaN());
                continue;
            }
            
            val = val.substr(start, end - start + 1);
            
            if (val == "NaN" || val == "nan" || val == "NAN") {
                raw_hits.push_back(std::numeric_limits<double>::quiet_NaN());
            } else {
                try {
                    raw_hits.push_back(std::stod(val));
                } catch (const std::exception &e) {
                    std::cerr << "ERROR: Cannot parse value '" << val << "' as double: " << e.what() << std::endl;
                    raw_hits.push_back(std::numeric_limits<double>::quiet_NaN());
                }
            }
        }
        if (raw_hits.size() < 8)
            continue;

        // --- HIT CORRECTION ---
        // Subtract the shift from the hits before the fit (only for active layers)
        std::vector<double> corrected_hits(8, 0.0);
        for (int i = 0; i < 8; ++i)
        {
            if (!IS_ACTIVE[i]) {
                corrected_hits[i] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }
            int shift_label = GLOBAL_LABELS[i][0];                 // shift ID
            double current_shift = alignment_offsets[shift_label];
            corrected_hits[i] = raw_hits[i] - current_shift;
        }

        // --- PRE-SELECTION: keep events in the main residual peak ---
        if(iteration == 1)
        {
            cY_ref = 0.085;   // main-peak position from Residual_0
            cX_ref = -0.04;    // read this off Residual_3 (or 5/7) for the X view
        }
        else
        {
            cY_ref = 0.015;   // main-peak position from Residual_0
            cX_ref = -0.01;    // read this off Residual_3 (or 5/7) for the X view
        }
        double dY = corrected_hits[0] - 2*corrected_hits[2] + corrected_hits[4];
        double dX = corrected_hits[3] - 2*corrected_hits[5] + corrected_hits[7];
        if (std::isnan(dX) || std::isnan(dY)) continue;
        if (std::abs(dY/6 - cY_ref) > CUT || std::abs(dX/6 - cX_ref) > CUT) continue;

        // --- PHASE 1: FIT LOCALLY WITH CORRECTED HITS ---
        fitterX.ClearPoints();
        fitterY.ClearPoints();

        for (int i : idxX)
        {
            if (!IS_ACTIVE[i] || std::isnan(corrected_hits[i])) continue;
            double z = Z_LAYERS[i];
            fitterX.AddPoint(&z, corrected_hits[i]);
        }
        for (int i : idxY)
        {
            if (!IS_ACTIVE[i] || std::isnan(corrected_hits[i])) continue;
            double z = Z_LAYERS[i];
            fitterY.AddPoint(&z, corrected_hits[i]);
        }

        fitterX.Eval();
        fitterY.Eval();

        double x0 = fitterX.GetParameter(0);
        double tx = fitterX.GetParameter(1);
        double y0 = fitterY.GetParameter(0);
        double ty = fitterY.GetParameter(1);

        // --- PHASE 2: CALCULATE THE RESIDUAL AND DERIVATIVES ---
        for (int i = 0; i < 8; ++i)
        {
            if (!IS_ACTIVE[i] || std::isnan(corrected_hits[i])) continue;

            double z = Z_LAYERS[i];
            double measurement = corrected_hits[i];
            double x_pred = x0 + tx * z;
            double y_pred = y0 + ty * z;

            float residual = 0.0;
            float derLc[4] = {0.0, 0.0, 0.0, 0.0};
            float derGl[4] = {0.0, 0.0, 0.0, 0.0};

            if (IS_X_LAYER[i])
            {
                residual = static_cast<float>(measurement - x_pred);
                derLc[0] = 1.0f;
                derLc[1] = static_cast<float>(z);
                derGl[0] = 1.0f;
                derGl[1] = static_cast<float>(-y_pred);
                derGl[2] = static_cast<float>(y_pred * tx);
                derGl[3] = static_cast<float>(-z * tx);
            }
            else
            {
                residual = static_cast<float>(measurement - y_pred);
                derLc[2] = 1.0f;
                derLc[3] = static_cast<float>(z);
                derGl[0] = 1.0f;
                derGl[1] = static_cast<float>(x_pred);
                derGl[2] = static_cast<float>(z * ty);
                derGl[3] = static_cast<float>(-x_pred * ty);
            }

            residuals[i]->Fill(residual);

            milleWriter.mille(4, derLc, 4, derGl, GLOBAL_LABELS[i], residual, SIGMA_HIT);
        }

        milleWriter.end();
        eventCount++;
    }

    outputFile.cd();
    for (int i = 0; i < 8; ++i)
    {
        residuals[i]->Write();
    }
    outputFile.Close();
    
    milleWriter.kill();
    return 0;
}