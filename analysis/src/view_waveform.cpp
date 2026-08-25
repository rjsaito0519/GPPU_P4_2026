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
#include <Buttons.h>

using namespace std;

static const int _DT5751Length = 1029;

// イベント制御用のグローバル変数
TCanvas* c_wave = nullptr;
TTree* wave_tree = nullptr;
std::vector<Int_t> event_ids;
std::map<Int_t, std::pair<Long64_t, Long64_t>> event_map;
size_t current_idx = 0;
TH1D* h_ch0 = nullptr;
TH1D* h_ch1 = nullptr;

Int_t g_event;
Int_t g_channel;
ULong64_t g_time_stamp;
UShort_t g_wave_raw[_DT5751Length];

// 指定されたインデックスのイベントを描画する関数
void DrawEvent(size_t idx) {
    if (idx >= event_ids.size()) return;
    Int_t ev_id = event_ids[idx];
    auto entries = event_map[ev_id];

    c_wave->cd(1);
    gPad->Clear();
    if (entries.first != -1) {
        wave_tree->GetEntry(entries.first);
        h_ch0->Reset();
        h_ch0->SetTitle(Form("CH0 - Event %d (timestamp: %llu)", ev_id, (unsigned long long)g_time_stamp));
        double ymin = 99999;
        double ymax = -99999;
        for (int k = 0; k < _DT5751Length; k++) {
            h_ch0->SetBinContent(k + 1, g_wave_raw[k]);
            if (g_wave_raw[k] < ymin) ymin = g_wave_raw[k];
            if (g_wave_raw[k] > ymax) ymax = g_wave_raw[k];
        }
        h_ch0->GetYaxis()->SetRangeUser(ymin - 20.0, ymax + 20.0);
        h_ch0->Draw("hist");
    } else {
        TH1D* h_empty = new TH1D("h_empty0", Form("CH0 - Event %d (No Data);time [ns];ADC Value", ev_id), _DT5751Length, 0, _DT5751Length);
        h_empty->Draw("hist");
    }

    c_wave->cd(2);
    gPad->Clear();
    if (entries.second != -1) {
        wave_tree->GetEntry(entries.second);
        h_ch1->Reset();
        h_ch1->SetTitle(Form("CH1 - Event %d (timestamp: %llu)", ev_id, (unsigned long long)g_time_stamp));
        double ymin = 99999;
        double ymax = -99999;
        for (int k = 0; k < _DT5751Length; k++) {
            h_ch1->SetBinContent(k + 1, g_wave_raw[k]);
            if (g_wave_raw[k] < ymin) ymin = g_wave_raw[k];
            if (g_wave_raw[k] > ymax) ymax = g_wave_raw[k];
        }
        h_ch1->GetYaxis()->SetRangeUser(ymin - 20.0, ymax + 20.0);
        h_ch1->Draw("hist");
    } else {
        TH1D* h_empty = new TH1D("h_empty1", Form("CH1 - Event %d (No Data);time [ns];ADC Value", ev_id), _DT5751Length, 0, _DT5751Length);
        h_empty->Draw("hist");
    }

    c_wave->Update();
    cout << Form("\r[Event ID: %d] index %d/%d (Press 'n': Next, 'p': Prev, 'q': Quit on Canvas)" , ev_id, (int)idx, (int)(event_ids.size() - 1)) << flush;
}

// キャンバス上でのキーボード入力を処理するハンドラ関数
void HandleKeyPress(Int_t event, Int_t x, Int_t y, TObject* selected) {
    // event == 24 は ROOT における kKeyPress イベント (X11上のキー押下)
    if (event == kKeyPress) {
        // x に押されたキーの ASCIIコードが入ります
        if (x == 'n' || x == 'N') {
            if (current_idx + 1 < event_ids.size()) {
                current_idx++;
                DrawEvent(current_idx);
            } else {
                cout << "\nAlready at the last event." << endl;
            }
        } else if (x == 'p' || x == 'P') {
            if (current_idx > 0) {
                current_idx--;
                DrawEvent(current_idx);
            } else {
                cout << "\nAlready at the first event." << endl;
            }
        } else if (x == 'q' || x == 'Q') {
            cout << "\nTerminating viewer..." << endl;
            gApplication->Terminate(0);
        }
    }
}

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

    wave_tree = (TTree*)file->Get("wave_tree");
    if (!wave_tree) {
        cerr << "ERROR: cannot find TTree 'wave_tree' in input file" << endl;
        file->Close();
        return 1;
    }

    // インタラクティブ描画のために TApplication を初期化
    TApplication app("app", &argc, argv);

    wave_tree->SetBranchAddress("event", &g_event);
    wave_tree->SetBranchAddress("channel", &g_channel);
    wave_tree->SetBranchAddress("time_stamp", &g_time_stamp);
    wave_tree->SetBranchAddress("wave_raw", g_wave_raw);

    // イベントマップの構築 (event -> {ch0_entry, ch1_entry})
    Long64_t n_entries = wave_tree->GetEntries();
    cout << "Scanning TTree and building event index..." << endl;
    for (Long64_t i = 0; i < n_entries; ++i) {
        wave_tree->GetEntry(i);
        if (event_map.find(g_event) == event_map.end()) {
            event_map[g_event] = make_pair(-1LL, -1LL);
            event_ids.push_back(g_event);
        }
        if (g_channel == 0) {
            event_map[g_event].first = i;
        } else if (g_channel == 1) {
            event_map[g_event].second = i;
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
    c_wave = new TCanvas("c_wave", "Waveform Viewer", 800, 800);
    c_wave->Divide(1, 2);

    // キャンバスのイベント信号（キー押下など）を HandleKeyPress グローバル関数に接続
    c_wave->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", nullptr, nullptr, "HandleKeyPress(Int_t,Int_t,Int_t,TObject*)");

    // 描画用ヒストグラムの作成
    h_ch0 = new TH1D("h_ch0", "CH0 Raw Waveform;time [ns];ADC Value", _DT5751Length, 0, _DT5751Length);
    h_ch1 = new TH1D("h_ch1", "CH1 Raw Waveform;time [ns];ADC Value", _DT5751Length, 0, _DT5751Length);
    h_ch0->SetLineColor(kBlue);
    h_ch1->SetLineColor(kRed);
    h_ch0->SetLineWidth(1);
    h_ch1->SetLineWidth(1);

    // 初期イベントの描画
    current_idx = 0;
    DrawEvent(current_idx);

    // ROOT GUI のメインイベントループを開始 (キー入力を監視して待機)
    app.Run();

    file->Close();
    delete file;
    return 0;
}
