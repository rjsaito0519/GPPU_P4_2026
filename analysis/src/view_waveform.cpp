#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TApplication.h>
#include <TSystem.h>

using namespace std;

static const Int_t _DT5751Length = 1029;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_waveform_root_file>" << endl;
        return 1;
    }

    string input_path = argv[1];

    TFile* file = TFile::Open(input_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "ERROR: cannot open input ROOT file -> " << input_path << endl;
        return 1;
    }

    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        cerr << "ERROR: cannot find TTree 'tree' in input file" << endl;
        file->Close();
        return 1;
    }

    // インタラクティブ描画のために TApplication を初期化
    TApplication app("app", &argc, argv);

    Int_t event;
    Int_t channel;
    ULong64_t time_stamp;
    UShort_t wave_raw[_DT5751Length];

    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("channel", &channel);
    tree->SetBranchAddress("time_stamp", &time_stamp);
    tree->SetBranchAddress("wave_raw", wave_raw);

    // イベントマップの構築 (event -> {ch0_entry, ch1_entry})
    map<Int_t, pair<Long64_t, Long64_t>> event_map;
    vector<Int_t> event_ids;

    Long64_t n_entries = tree->GetEntries();
    cout << "Scanning TTree and building event index..." << endl;
    for (Long64_t i = 0; i < n_entries; ++i) {
        tree->GetEntry(i);
        if (event_map.find(event) == event_map.end()) {
            event_map[event] = make_pair(-1LL, -1LL);
            event_ids.push_back(event);
        }
        if (channel == 0) {
            event_map[event].first = i;
        } else if (channel == 1) {
            event_map[event].second = i;
        }
    }

    // イベントIDをソート
    sort(event_ids.begin(), event_ids.end());

    if (event_ids.empty()) {
        cerr << "ERROR: No events found in tree." << endl;
        file->Close();
        return 1;
    }

    cout << "Index build complete. Unique Events: " << event_ids.size() << endl;

    // キャンバスのセットアップ (上下2画面)
    TCanvas* c_wave = new TCanvas("c_wave", "Waveform Viewer", 800, 800);
    c_wave->Divide(1, 2);

    // 描画用ヒストグラムの作成
    TH1D* h_ch0 = new TH1D("h_ch0", "CH0 Raw Waveform;time [ns];ADC Value", _DT5751Length, 0, _DT5751Length);
    TH1D* h_ch1 = new TH1D("h_ch1", "CH1 Raw Waveform;time [ns];ADC Value", _DT5751Length, 0, _DT5751Length);
    h_ch0->SetLineColor(kBlue);
    h_ch1->SetLineColor(kRed);
    h_ch0->SetLineWidth(1);
    h_ch1->SetLineWidth(1);

    size_t current_idx = 0;

    auto DrawEvent = [&](size_t idx) {
        if (idx >= event_ids.size()) return;
        Int_t ev_id = event_ids[idx];
        auto entries = event_map[ev_id];

        c_wave->cd(1);
        gPad->Clear();
        if (entries.first != -1) {
            tree->GetEntry(entries.first);
            h_ch0->Reset();
            h_ch0->SetTitle(Form("CH0 - Event %d (timestamp: %llu)", ev_id, (unsigned long long)time_stamp));
            double ymin = 99999;
            double ymax = -99999;
            for (int k = 0; k < _DT5751Length; k++) {
                h_ch0->SetBinContent(k + 1, wave_raw[k]);
                if (wave_raw[k] < ymin) ymin = wave_raw[k];
                if (wave_raw[k] > ymax) ymax = wave_raw[k];
            }
            h_ch0->GetYaxis()->SetRangeUser(ymin - 20.0, ymax + 20.0);
            h_ch0->Draw("hist");
        } else {
            // CH0データがない場合
            TH1D* h_empty = new TH1D("h_empty0", Form("CH0 - Event %d (No Data);time [ns];ADC Value", ev_id), _DT5751Length, 0, _DT5751Length);
            h_empty->Draw("hist");
        }

        c_wave->cd(2);
        gPad->Clear();
        if (entries.second != -1) {
            tree->GetEntry(entries.second);
            h_ch1->Reset();
            h_ch1->SetTitle(Form("CH1 - Event %d (timestamp: %llu)", ev_id, (unsigned long long)time_stamp));
            double ymin = 99999;
            double ymax = -99999;
            for (int k = 0; k < _DT5751Length; k++) {
                h_ch1->SetBinContent(k + 1, wave_raw[k]);
                if (wave_raw[k] < ymin) ymin = wave_raw[k];
                if (wave_raw[k] > ymax) ymax = wave_raw[k];
            }
            h_ch1->GetYaxis()->SetRangeUser(ymin - 20.0, ymax + 20.0);
            h_ch1->Draw("hist");
        } else {
            // CH1データがない場合
            TH1D* h_empty = new TH1D("h_empty1", Form("CH1 - Event %d (No Data);time [ns];ADC Value", ev_id), _DT5751Length, 0, _DT5751Length);
            h_empty->Draw("hist");
        }

        c_wave->Update();
        gSystem->ProcessEvents(); // X11のグラフィックバッファを強制フラッシュ
    };

    // 初期イベントの描画
    DrawEvent(current_idx);

    // インタラクティブ入力ループ
    string cmd;
    cout << "\n========================================================" << endl;
    cout << " Waveform Viewer Interactive Mode (Console Control)" << endl;
    cout << "   n + Enter : Go to Next Event" << endl;
    cout << "   p + Enter : Go to Previous Event" << endl;
    cout << "   [Event ID] + Enter : Jump directly to specific Event ID" << endl;
    cout << "   q + Enter : Quit Viewer" << endl;
    cout << "========================================================" << endl;

    while (true) {
        Int_t ev_id = event_ids[current_idx];
        cout << Form("\n[Current Event ID: %d] index %d/%d >> ", ev_id, (int)current_idx, (int)(event_ids.size() - 1)) << flush;
        if (!(cin >> cmd)) break;

        if (cmd == "q" || cmd == "Q") {
            break;
        } else if (cmd == "n" || cmd == "N") {
            if (current_idx + 1 < event_ids.size()) {
                current_idx++;
                DrawEvent(current_idx);
            } else {
                cout << "Already at the last event." << endl;
            }
        } else if (cmd == "p" || cmd == "P") {
            if (current_idx > 0) {
                current_idx--;
                DrawEvent(current_idx);
            } else {
                cout << "Already at the first event." << endl;
            }
        } else {
            // 直接のイベントIDジャンプを試行
            try {
                int jump_ev = stoi(cmd);
                auto it = find(event_ids.begin(), event_ids.end(), jump_ev);
                if (it != event_ids.end()) {
                    current_idx = distance(event_ids.begin(), it);
                    DrawEvent(current_idx);
                } else {
                    cout << "Event ID " << jump_ev << " not found in this file." << endl;
                }
            } catch (...) {
                cout << "Invalid command. Use 'n', 'p', 'q', or numeric Event ID." << endl;
            }
        }
    }

    file->Close();
    delete file;
    return 0;
}
