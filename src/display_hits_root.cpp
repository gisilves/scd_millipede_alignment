#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TStyle.h"

constexpr int N_DETECTORS = 8;
constexpr float DETECTOR_WIDTH = 200.0f;

struct HitData {
    float x, y, z;
    int detector_id;
    bool is_x_layer;
};

std::vector<HitData> readHits(const std::string& csv_file) {
    std::vector<HitData> hits;
    
    std::ifstream csv(csv_file);
    if (!csv.is_open()) {
        std::cerr << "ERROR: Cannot open " << csv_file << std::endl;
        return hits;
    }
    
    std::string line;
    int event_id = 0;
    
    while (std::getline(csv, line)) {
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string token;
        int det = 0;
        
        while (std::getline(ss, token, ',') && det < N_DETECTORS) {
            if (!token.empty() && token != "NaN" && token != "nan") {
                try {
                    float val = std::stof(token);
                    
                    // Reconstruct position
                    float z = det * 10.0f;  // Layer spacing
                    bool is_x = (det % 2 == 0);
                    
                    HitData hit;
                    hit.z = z;
                    hit.detector_id = det;
                    hit.is_x_layer = is_x;
                    
                    if (is_x) {
                        hit.x = val;
                        hit.y = 0.0f;  // Y not measured on X-layer
                    } else {
                        hit.x = 0.0f;  // X not measured on Y-layer
                        hit.y = val;
                    }
                    
                    hits.push_back(hit);
                    
                } catch (...) {
                    // Skip invalid entries
                }
            }
            det++;
        }
        
        event_id++;
    }
    
    csv.close();
    return hits;
}

void createHistograms(const std::vector<HitData>& hits, 
                     const std::string& output_file) {
    
    // Create ROOT file
    TFile* file = TFile::Open(output_file.c_str(), "RECREATE");
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: Cannot create " << output_file << std::endl;
        return;
    }
    
    std::cout << "Creating histograms with " << hits.size() << " hits..." << std::endl;
    
    // 3D histogram: x vs y vs z (detector position space)
    TH3F* h3_xyz = new TH3F("h3_xyz", 
                            "Hit Distribution (XYZ);X (mm);Y (mm);Z (mm)",
                            100, -20, 20, 100, -20, 20, 8, -5, 75);
    h3_xyz->Sumw2();
    
    // 2D histogram: x vs z
    TH2F* h2_xz = new TH2F("h2_xz",
                           "Hit Distribution (XZ);X (mm);Z (mm)",
                           50, -20, 20, 8, -5, 75);
    h2_xz->Sumw2();
    
    // 2D histogram: y vs z
    TH2F* h2_yz = new TH2F("h2_yz",
                           "Hit Distribution (YZ);Y (mm);Z (mm)",
                           50, -20, 20, 8, -5, 75);
    h2_yz->Sumw2();
    
    // Per-detector histograms
    std::array<TH1F*, N_DETECTORS> h1_det;
    for (int det = 0; det < N_DETECTORS; det++) {
        std::string name = std::string("h1_det_") + std::to_string(det);
        std::string title = std::string("Detector ") + std::to_string(det) + 
                           (det % 2 == 0 ? " (X-layer)" : " (Y-layer)");
        h1_det[det] = new TH1F(name.c_str(), title.c_str(), 100, -20, 20);
        h1_det[det]->Sumw2();
        h1_det[det]->SetXTitle("Position (mm)");
        h1_det[det]->SetYTitle("Count");
    }
    
    // TTree for detailed event-by-event data
    TTree* tree = new TTree("hits", "Hit data tree");
    
    float x, y, z;
    int det_id;
    bool is_x;
    
    tree->Branch("x", &x, "x/F");
    tree->Branch("y", &y, "y/F");
    tree->Branch("z", &z, "z/F");
    tree->Branch("detector_id", &det_id, "detector_id/I");
    tree->Branch("is_x_layer", &is_x, "is_x_layer/O");
    
    // Fill histograms
    int n_x_hits = 0, n_y_hits = 0;
    
    for (const auto& hit : hits) {
        // 3D histogram
        if (hit.is_x_layer) {
            h3_xyz->Fill(hit.x, 0.0, hit.z);
            h2_xz->Fill(hit.x, hit.z);
            n_x_hits++;
        } else {
            h3_xyz->Fill(0.0, hit.y, hit.z);
            h2_yz->Fill(hit.y, hit.z);
            n_y_hits++;
        }
        
        // Per-detector histogram
        float pos = hit.is_x_layer ? hit.x : hit.y;
        h1_det[hit.detector_id]->Fill(pos);
        
        // TTree
        x = hit.x;
        y = hit.y;
        z = hit.z;
        det_id = hit.detector_id;
        is_x = hit.is_x_layer;
        tree->Fill();
    }
    
    std::cout << "Filled: " << n_x_hits << " X-hits, " << n_y_hits 
              << " Y-hits" << std::endl;
    
    // Write histograms
    file->cd();
    h3_xyz->Write();
    h2_xz->Write();
    h2_yz->Write();
    tree->Write();
    
    for (int det = 0; det < N_DETECTORS; det++) {
        h1_det[det]->Write();
    }
    
    file->Close();
    std::cout << "Wrote histograms to " << output_file << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <hits.csv> [output.root]" << std::endl;
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string output_file = (argc > 2) ? argv[2] : "hits_histograms.root";
    
    std::cout << "Reading hits from " << input_file << std::endl;
    auto hits = readHits(input_file);
    
    if (hits.empty()) {
        std::cerr << "ERROR: No hits found" << std::endl;
        return 1;
    }
    
    createHistograms(hits, output_file);
    
    return 0;
}
