#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include "TLinearFitter.h"
#include "Mille.h"
#include "TFile.h"
#include "TH1.h"

const std::vector<double> Z_LAYERS = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0};
const std::vector<bool> IS_X_LAYER = {true, false, true, false, true, false, true, false};
const float SIGMA_HIT = 0.05;

const int GLOBAL_LABELS[8][4] ={
    {1011, 1013, 1014, 1015}, // Strato 1 (X)
    {2012, 2013, 2014, 2015}, // Strato 2 (Y)
    {3011, 3013, 3014, 3015}, // Strato 3 (X)
    {4012, 4013, 4014, 4015}, // Strato 4 (Y)
    {5011, 5013, 5014, 5015}, // Strato 5 (X)
    {6012, 6013, 6014, 6015}, // Strato 6 (Y)
    {7011, 7013, 7014, 7015}, // Strato 7 (X)
    {8012, 8013, 8014, 8015}  // Strato 8 (Y)
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

int main()
{
    // 1. Load the alignment corrections calculated in previous steps
    loadGeometry();

    std::ifstream csvFile("hits_generated.csv");
    if (!csvFile.is_open())
    {
        std::cerr << "Error: cannot open the CSV file with the measurements." << std::endl;
        return 1;
    }

    // 2. Output file for residual histograms
    TFile outputFile("residuals.root", "RECREATE");
    std::vector<TH1*> residuals;
    for (int i = 0; i < 8; ++i)
    {
        std::string name = "Residual_" + std::to_string(i);
        residuals.push_back(new TH1F(name.c_str(), name.c_str(), 500, -1.0, 1.0));
    }

    Mille milleWriter("mille_output_run001.bin");
    std::string line;
    int eventCount = 0;

    TLinearFitter fitterX(1, "pol1");
    TLinearFitter fitterY(1, "pol1");

    std::vector<int> idxX = {0, 2, 4, 6};
    std::vector<int> idxY = {1, 3, 5, 7};

    while (std::getline(csvFile, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::string val;
        std::vector<double> raw_hits;

        while (std::getline(ss, val, ','))
        {
            raw_hits.push_back(std::stod(val));
        }
        if (raw_hits.size() < 8)
            continue;

        // --- HIT CORRECTION ---
        // Subtract the shift (xx11 or xx12) from the hits before the fit
        std::vector<double> corrected_hits(8, 0.0);
        for (int i = 0; i < 8; ++i)
        {
            int shift_label = GLOBAL_LABELS[i][0];                 // shift ID (eg. 3011 or 4012)
            double current_shift = alignment_offsets[shift_label];
            corrected_hits[i] = raw_hits[i] - current_shift;
        }

        // --- PHASE 1: FIT LOCALLY WITH CORRECTED HITS ---
        fitterX.ClearPoints();
        fitterY.ClearPoints();

        for (int i : idxX)
        {
            double z = Z_LAYERS[i];
            fitterX.AddPoint(&z, corrected_hits[i]);
        }
        for (int i : idxY)
        {
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
