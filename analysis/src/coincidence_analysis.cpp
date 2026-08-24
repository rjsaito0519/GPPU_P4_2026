#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include "progress_bar.h"

// メモリ上にロードするためのイベントデータ構造体
struct EventData {
    Int_t event;
    Long64_t TS;
    Double_t T0;
    Double_t Q0;
    Double_t T1;
    Double_t Q1;
    Double_t T0_corr;
    Int_t data_id;
};

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_root_path> [output_root_path]" << std::endl;
        return 1;
    }

    std::string input_path = argv[1];
    std::string output_path;
    if (argc > 2) {
        output_path = argv[2];
    } else {
        size_t last_dot = input_path.find_last_of(".");
        if (last_dot != std::string::npos) {
            output_path = input_path.substr(0, last_dot) + "_coincidence.root";
        } else {
            output_path = input_path + "_coincidence.root";
        }
    }

    TFile* fin = TFile::Open(input_path.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << input_path << std::endl;
        return 1;
    }

    // TTreeReader でデータをメモリ上にロードする
    TTreeReader reader("tree", fin);
    TTreeReaderValue<Int_t> event(reader, "event");
    TTreeReaderValue<Long64_t> TS(reader, "TS");
    TTreeReaderValue<Double_t> T0(reader, "T0");
    TTreeReaderValue<Double_t> Q0(reader, "Q0");
    TTreeReaderValue<Double_t> T1(reader, "T1");
    TTreeReaderValue<Double_t> Q1(reader, "Q1");
    TTreeReaderValue<Double_t> T0_corr(reader, "T0_corr");
    TTreeReaderValue<Int_t> data_id(reader, "data_id");

    Long64_t nentries = reader.GetEntries(true);
    std::cout << "Loading " << nentries << " entries into memory..." << std::endl;

    std::vector<EventData> events;
    events.reserve(nentries);

    Long64_t load_count = 0;
    while (reader.Next()) {
        if (load_count % 50000 == 0 || load_count == nentries - 1) {
            displayProgressBar(load_count + 1, nentries);
        }
        events.push_back({*event, *TS, *T0, *Q0, *T1, *Q1, *T0_corr, *data_id});
        load_count++;
    }
    std::cout << "\nLoaded " << events.size() << " entries successfully." << std::endl;

    // 出力ファイルの準備
    TFile* fout = new TFile(output_path.c_str(), "RECREATE");
    TTree* out_tree = new TTree("coincidence_tree", "Coincidence Analysis Results");

    Int_t out_fast_event;
    Long64_t out_fast_TS;
    Int_t out_slow_event;
    Long64_t out_slow_TS;
    Long64_t out_delta_TS;
    Double_t out_delta_T_us;
    Double_t out_fast_T1;
    Double_t out_slow_T1;
    Double_t out_fast_Q1;
    Double_t out_slow_Q1;
    Double_t out_fast_T0_corr;
    Double_t out_slow_T0_corr;
    Int_t out_data_id;

    out_tree->Branch("fast_event", &out_fast_event, "fast_event/I");
    out_tree->Branch("fast_TS", &out_fast_TS, "fast_TS/L");
    out_tree->Branch("slow_event", &out_slow_event, "slow_event/I");
    out_tree->Branch("slow_TS", &out_slow_TS, "slow_TS/L");
    out_tree->Branch("delta_TS", &out_delta_TS, "delta_TS/L");
    out_tree->Branch("delta_T_us", &out_delta_T_us, "delta_T_us/D");
    out_tree->Branch("fast_T1", &out_fast_T1, "fast_T1/D");
    out_tree->Branch("slow_T1", &out_slow_T1, "slow_T1/D");
    out_tree->Branch("fast_Q1", &out_fast_Q1, "fast_Q1/D");
    out_tree->Branch("slow_Q1", &out_slow_Q1, "slow_Q1/D");
    out_tree->Branch("fast_T0_corr", &out_fast_T0_corr, "fast_T0_corr/D");
    out_tree->Branch("slow_T0_corr", &out_slow_T0_corr, "slow_T0_corr/D");
    out_tree->Branch("data_id", &out_data_id, "data_id/I");

    std::cout << "Analyzing coincidence events (using T0_corr for timing)..." << std::endl;
    Long64_t total_events = events.size();

    for (Long64_t i = 0; i < total_events; ++i) {
        if (i % 5000 == 0 || i == total_events - 1) {
            displayProgressBar(i + 1, total_events);
        }

        const auto& fast = events[i];
        // 補正された T0_corr を用いて TOF を判定する
        Double_t t_diff = fast.T1 - fast.T0_corr;

        // 起点(fast)の条件チェック: T1-T0_corr が 60 ~ 120 ns の間
        if (t_diff >= 60.0 && t_diff <= 120.0) {
            // 後続イベントを1 ms (125000 TS) 以内でスキャン
            for (Long64_t j = i + 1; j < total_events; ++j) {
                const auto& slow = events[j];

                // 異なる data_id に到達したら、このファイルのデータは終了なので探索を打ち切る
                if (slow.data_id != fast.data_id) {
                    break;
                }

                Long64_t diff_TS = slow.TS - fast.TS;

                // 125000 TS を超えたら探索終了
                if (diff_TS > 125000) {
                    break;
                }

                // 後発(slow)のイベントで T1 ヒットがある（T1 > 0）ものをすべて保存
                if (slow.T1 > 0.0) {
                    out_fast_event = fast.event;
                    out_fast_TS = fast.TS;
                    out_slow_event = slow.event;
                    out_slow_TS = slow.TS;
                    out_delta_TS = diff_TS;
                    out_delta_T_us = diff_TS * 0.008; // 1 TS = 8 ns = 0.008 us
                    out_fast_T1 = fast.T1;
                    out_slow_T1 = slow.T1;
                    out_fast_Q1 = fast.Q1;
                    out_slow_Q1 = slow.Q1;
                    out_fast_T0_corr = fast.T0_corr;
                    out_slow_T0_corr = slow.T0_corr;
                    out_data_id = fast.data_id;

                    out_tree->Fill();
                }
            }
        }
    }

    std::cout << "\nAnalysis completed. Results saved to " << output_path << std::endl;

    out_tree->Write();
    fout->Close();
    fin->Close();

    delete fout;
    delete fin;

    return 0;
}
