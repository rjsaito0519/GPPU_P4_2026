#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLegend.h>

using namespace std;

static const int _DT5751Length = 1029;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_root_file> [--n num_events] [--ch channel] [--norm 0|1]" << endl;
        return 1;
    }

    string input_path = argv[1];
    int num_to_overlay = 20;
    int target_ch = 1;
    int norm_flag = 1; // デフォルトで縦軸を 1.0 に規格化する

    for (int i = 2; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--n" && i + 1 < argc) {
            num_to_overlay = stoi(argv[++i]);
        } else if (arg == "--ch" && i + 1 < argc) {
            target_ch = stoi(argv[++i]);
        } else if (arg == "--norm" && i + 1 < argc) {
            norm_flag = stoi(argv[++i]);
        }
    }

    TFile* file = TFile::Open(input_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "ERROR: cannot open input ROOT file -> " << input_path << endl;
        return 1;
    }

    TTree* wave_tree = (TTree*)file->Get("wave_tree");
    if (!wave_tree) {
        cerr << "ERROR: cannot find TTree 'wave_tree' in input file" << endl;
        file->Close();
        return 1;
    }

    // インタラクティブ画面表示用の TApplication
    TApplication app("app", &argc, argv);

    Int_t event;
    Int_t channel;
    ULong64_t time_stamp;
    UShort_t wave_raw[_DT5751Length];

    wave_tree->SetBranchAddress("event", &event);
    wave_tree->SetBranchAddress("channel", &channel);
    wave_tree->SetBranchAddress("time_stamp", &time_stamp);
    wave_tree->SetBranchAddress("wave_raw", wave_raw);

    TCanvas* c = new TCanvas("c_overlay", "Waveform Overlay Plot", 900, 700);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetGrid();

    TMultiGraph* mg = new TMultiGraph();
    string title = "Superimposed Waveforms (Aligned at Peak)";
    if (norm_flag) title += " [Normalized]";
    title += ";Time relative to peak [ns];" + string(norm_flag ? "Normalized Amplitude" : "Pulse Height [ADC]");
    mg->SetTitle(title.c_str());

    Long64_t n_entries = wave_tree->GetEntries();
    int count = 0;

    const double threshold = 10.0;
    const int n_baseline_length = 200;

    // カラーパレットの準備 (重ね書き用に見やすい色を選択)
    vector<int> colors = {kBlue, kRed, kGreen+2, kOrange+1, kMagenta, kCyan+2, kViolet, kPink-3, kTeal+9, kAzure+2};

    cout << "Scanning events and aligning waveforms..." << endl;

    for (Long64_t i = 0; i < n_entries; ++i) {
        if (count >= num_to_overlay) break;

        wave_tree->GetEntry(i);

        if (channel != target_ch) {
            continue;
        }

        // 1. ベースライン計算
        double sum_wave = 0.0;
        int n_wave = 0;
        for (int k = 0; k < n_baseline_length; ++k) {
            sum_wave += (double)wave_raw[k];
            n_wave++;
        }
        double baseline_rough = sum_wave / double(n_wave);

        sum_wave = 0.0;
        n_wave = 0;
        for (int k = 0; k < n_baseline_length; ++k) {
            if (abs((double)wave_raw[k] - baseline_rough) < threshold) {
                sum_wave += (double)wave_raw[k];
                n_wave++;
            }
        }
        double baseline = (n_wave > 0) ? (sum_wave / double(n_wave)) : baseline_rough;

        // 2. 反転差分波形の作成
        vector<double> wave(_DT5751Length);
        for (int k = 0; k < _DT5751Length; ++k) {
            wave[k] = baseline - (double)wave_raw[k];
        }

        // 3. ピークサーチ (300 ~ 500 ns)
        double max_val = -99999.0;
        int k_peak = -1;
        for (int k = 300; k < 500; ++k) {
            if (wave[k] > max_val) {
                max_val = wave[k];
                k_peak = k;
            }
        }

        if (k_peak == -1 || max_val < threshold) {
            continue;
        }

        // 4. ピーク合わせ TGraph の構築
        vector<double> x(_DT5751Length);
        vector<double> y(_DT5751Length);

        for (int k = 0; k < _DT5751Length; ++k) {
            x[k] = (double)(k - k_peak); // 横軸: ピーク位置を 0 ns に揃える
            if (norm_flag) {
                y[k] = wave[k] / max_val; // 縦軸: ピーク高を 1.0 に規格化
            } else {
                y[k] = wave[k]; // 縦軸: 生のパルス高 (ADC)
            }
        }

        TGraph* g = new TGraph(_DT5751Length, &x[0], &y[0]);
        int col = colors[count % colors.size()];
        g->SetLineColor(col);
        g->SetLineWidth(1);
        
        mg->Add(g, "L");
        count++;
    }

    if (count == 0) {
        cerr << "WARNING: No valid waveforms found matching channel " << target_ch << endl;
        file->Close();
        return 1;
    }

    cout << "Overlayed " << count << " waveforms." << endl;

    mg->Draw("A");
    
    // 表示レンジの調整 (ピーク前後を表示)
    mg->GetXaxis()->SetRangeUser(-20.0, 150.0);
    if (norm_flag) {
        mg->GetYaxis()->SetRangeUser(-0.1, 1.1);
    }
    
    c->Update();
    
    // 画像およびPDFとして自動保存
    string out_pdf = "pdf/waveform_overlay_ch" + to_string(target_ch) + (norm_flag ? "_norm" : "") + ".pdf";
    gSystem->mkdir("pdf", true);
    c->Print(out_pdf.c_str());
    cout << "Overlay plot saved to: " << out_pdf << endl;

    // ウインドウ表示を維持
    app.Run();

    file->Close();
    delete file;
    return 0;
}
