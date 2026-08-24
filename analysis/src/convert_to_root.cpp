#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include "progress_bar.h"

Long64_t count_lines(const std::string& filename) {
    std::ifstream ifs(filename);
    Long64_t lines = 0;
    std::string line;
    while (std::getline(ifs, line)) {
        lines++;
    }
    return lines;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_dat_path> <data_id> <p0_slew> [output_root_path]" << std::endl;
        return 1;
    }

    std::string input_path = argv[1];
    Int_t input_data_id = 1;
    Double_t p0_slew = 0.0;
    std::string output_path;

    try {
        input_data_id = std::stoi(argv[2]);
    } catch (const std::exception& e) {
        std::cerr << "Error: Invalid data_id provided." << std::endl;
        return 1;
    }

    try {
        p0_slew = std::stod(argv[3]);
    } catch (const std::exception& e) {
        std::cerr << "Error: Invalid p0_slew provided." << std::endl;
        return 1;
    }

    if (argc > 4) {
        output_path = argv[4];
    } else {
        // 入力ファイルの拡張子を.rootに変更
        size_t last_dot = input_path.find_last_of(".");
        if (last_dot != std::string::npos) {
            output_path = input_path.substr(0, last_dot) + ".root";
        } else {
            output_path = input_path + ".root";
        }
    }

    Long64_t total_entries = count_lines(input_path);
    if (total_entries == 0) {
        std::cerr << "Error: Input file is empty or does not exist: " << input_path << std::endl;
        return 1;
    }

    std::ifstream ifs(input_path);
    if (!ifs.is_open()) {
        std::cerr << "Error: Cannot open input file " << input_path << std::endl;
        return 1;
    }

    TFile* fout = new TFile(output_path.c_str(), "RECREATE");
    TTree* tree = new TTree("tree", "Tree from dat file");

    Int_t event;
    Long64_t TS;
    Double_t T0;
    Double_t Q0;
    Double_t T1;
    Double_t Q1;
    Double_t T0_corr; // 補正後のT0
    Int_t out_data_id = input_data_id;

    // TTreeのブランチ設定
    tree->Branch("event", &event, "event/I");
    tree->Branch("TS", &TS, "TS/L");
    tree->Branch("T0", &T0, "T0/D");
    tree->Branch("Q0", &Q0, "Q0/D");
    tree->Branch("T1", &T1, "T1/D");
    tree->Branch("Q1", &Q1, "Q1/D");
    tree->Branch("T0_corr", &T0_corr, "T0_corr/D");
    tree->Branch("data_id", &out_data_id, "data_id/I");

    std::cout << "Converting " << input_path << " to " << output_path << std::endl;
    std::cout << " - data_id : " << out_data_id << std::endl;
    std::cout << " - p0_slew : " << p0_slew << " (Formula: T0_corr = T0 - p0/sqrt(Q0))" << std::endl;

    Long64_t count = 0;
    while (ifs >> event >> TS >> T0 >> Q0 >> T1 >> Q1) {
        // スルーイング補正の適用
        if (p0_slew != 0.0 && Q0 > 0.0) {
            T0_corr = T0 - (p0_slew / std::sqrt(Q0));
        } else {
            T0_corr = T0;
        }

        tree->Fill();
        count++;
        if (count % 5000 == 0 || count == total_entries) {
            displayProgressBar(count, total_entries);
        }
    }

    if (count < total_entries) {
        std::cerr << "\n[Warning] Stopped reading before the end of the file." << std::endl;
        std::cerr << "Expected: " << total_entries << " entries, but only loaded: " << count << std::endl;
        if (ifs.eof()) {
            std::cerr << "Reason: Reached EOF (End of File) unexpectedly." << std::endl;
        } else if (ifs.fail()) {
            std::cerr << "Reason: Data format mismatch or parse error at line " << count + 1 << "." << std::endl;
            ifs.clear();
            std::string failed_line;
            if (std::getline(ifs, failed_line)) {
                std::cerr << "Offending line content: \"" << failed_line << "\"" << std::endl;
            }
        }
    } else {
        std::cout << "\nTotal " << count << " entries loaded successfully." << std::endl;
    }

    tree->Write();
    fout->Close();

    delete fout;
    return 0;
}
