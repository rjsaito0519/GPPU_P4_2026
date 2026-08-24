#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include "progress_bar.h"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <output_root_path> <input_root_path1> <input_root_path2> ..." << std::endl;
        return 1;
    }

    std::string output_path = argv[1];
    std::vector<std::string> input_paths;
    for (int i = 2; i < argc; ++i) {
        input_paths.push_back(argv[i]);
    }

    TFile* fout = new TFile(output_path.c_str(), "RECREATE");
    TTree* out_tree = new TTree("tree", "Merged Tree");
    out_tree->SetDirectory(fout); // 出力ファイルを明示的に関連付け

    Int_t out_event;
    Long64_t out_TS;
    Double_t out_T0;
    Double_t out_Q0;
    Double_t out_T1;
    Double_t out_Q1;
    Double_t out_T0_corr; // マージ対象に追加
    Int_t out_data_id;    // マージ元の識別ID

    out_tree->Branch("event", &out_event, "event/I");
    out_tree->Branch("TS", &out_TS, "TS/L");
    out_tree->Branch("T0", &out_T0, "T0/D");
    out_tree->Branch("Q0", &out_Q0, "Q0/D");
    out_tree->Branch("T1", &out_T1, "T1/D");
    out_tree->Branch("Q1", &out_Q1, "Q1/D");
    out_tree->Branch("T0_corr", &out_T0_corr, "T0_corr/D");
    out_tree->Branch("data_id", &out_data_id, "data_id/I");

    for (size_t file_idx = 0; file_idx < input_paths.size(); ++file_idx) {
        std::string input_path = input_paths[file_idx];
        std::cout << "\nMerging file [" << file_idx + 1 << "/" << input_paths.size() << "]: " << input_path << std::endl;

        TFile* fin = TFile::Open(input_path.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            std::cerr << "Warning: Cannot open input file " << input_path << ". Skipping." << std::endl;
            if (fin) delete fin;
            continue;
        }

        TTreeReader reader("tree", fin);
        TTreeReaderValue<Int_t> event(reader, "event");
        TTreeReaderValue<Long64_t> TS(reader, "TS");
        TTreeReaderValue<Double_t> T0(reader, "T0");
        TTreeReaderValue<Double_t> Q0(reader, "Q0");
        TTreeReaderValue<Double_t> T1(reader, "T1");
        TTreeReaderValue<Double_t> Q1(reader, "Q1");
        TTreeReaderValue<Double_t> T0_corr(reader, "T0_corr");

        Long64_t nentries = reader.GetEntries(true);
        Long64_t count = 0;

        // ファイルのインデックス（1, 2, 3...）を data_id として割り振る
        out_data_id = static_cast<Int_t>(file_idx + 1);

        while (reader.Next()) {
            if (count % 50000 == 0 || count == nentries - 1) {
                displayProgressBar(count + 1, nentries);
            }

            out_event = *event;
            out_TS = *TS;
            out_T0 = *T0;
            out_Q0 = *Q0;
            out_T1 = *T1;
            out_Q1 = *Q1;
            out_T0_corr = *T0_corr;

            out_tree->Fill();
            count++;
        }

        fin->Close();
        delete fin;
    }

    std::cout << "\nAll files merged. Writing output to " << output_path << std::endl;
    fout->cd();
    out_tree->Write();
    fout->Close();
    delete fout;

    return 0;
}
