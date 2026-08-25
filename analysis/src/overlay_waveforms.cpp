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
#include <TAxis.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>

using namespace std;

static const Int_t _DT5751Length = 1029;
 
int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_root_file> [--n num_events] [--ch channel] [--norm 0|1] [--pre pre_ns] [--short short_ns] [--long long_ns]" << endl;
        return 1;
    }
 
    string input_path = argv[1];
    Int_t num_to_overlay = 20;
    Int_t target_ch = 1;
    Int_t norm_flag = 1; // デフォルトで縦軸を 1.0 に規格化する
    Int_t n_pre_peak = 10;
    Int_t n_post_peak_short = 10;
    Int_t n_post_peak_long = 30;
 
    for (Int_t i = 2; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--n" && i + 1 < argc) {
            num_to_overlay = stoi(argv[++i]);
        } else if (arg == "--ch" && i + 1 < argc) {
            target_ch = stoi(argv[++i]);
        } else if (arg == "--norm" && i + 1 < argc) {
            norm_flag = stoi(argv[++i]);
        } else if (arg == "--pre" && i + 1 < argc) {
            n_pre_peak = stoi(argv[++i]);
        } else if (arg == "--short" && i + 1 < argc) {
            n_post_peak_short = stoi(argv[++i]);
        } else if (arg == "--long" && i + 1 < argc) {
            n_post_peak_long = stoi(argv[++i]);
        }
    }

    TFile* file = TFile::Open(input_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "ERROR: cannot open input ROOT file -> " << input_path << endl;
        return 1;
    }

    TTree* wave_tree = (TTree*)file->Get("tree");
    if (!wave_tree) {
        cerr << "ERROR: cannot find TTree 'tree' in input file" << endl;
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
    Int_t count = 0;
 
    const Double_t threshold = 10.0;
    const Int_t n_baseline_length = 200;
 
    // カラーパレットの準備 (重ね書き用に見やすい色を選択)
    vector<Int_t> colors = {kBlue, kRed, kGreen+2, kOrange+1, kMagenta, kCyan+2, kViolet, kPink-3, kTeal+9, kAzure+2};

    cout << "Scanning events and aligning waveforms..." << endl;

    for (Long64_t i = 0; i < n_entries; ++i) {
        if (count >= num_to_overlay) break;
 
        wave_tree->GetEntry(i);

        if (channel != target_ch) {
            continue;
        }

        // 1. ベースライン計算
        Double_t sum_wave = 0.0;
        Int_t n_wave = 0;
        for (Int_t k = 0; k < n_baseline_length; ++k) {
            sum_wave += (Double_t)wave_raw[k];
            n_wave++;
        }
        Double_t baseline_rough = sum_wave / Double_t(n_wave);
 
        sum_wave = 0.0;
        n_wave = 0;
        for (Int_t k = 0; k < n_baseline_length; ++k) {
            if (abs((Double_t)wave_raw[k] - baseline_rough) < threshold) {
                sum_wave += (Double_t)wave_raw[k];
                n_wave++;
            }
        }
        Double_t baseline = (n_wave > 0) ? (sum_wave / Double_t(n_wave)) : baseline_rough;
 
        // 2. 反転差分波形の作成
        vector<Double_t> wave(_DT5751Length);
        for (Int_t k = 0; k < _DT5751Length; ++k) {
            wave[k] = baseline - (Double_t)wave_raw[k];
        }
 
        // 3. ピークサーチ (300 ~ 500 ns)
        Double_t max_val = -99999.0;
        Int_t k_peak = -1;
        for (Int_t k = 300; k < 500; ++k) {
            if (wave[k] > max_val) {
                max_val = wave[k];
                k_peak = k;
            }
        }

        if (k_peak == -1 || max_val < threshold) {
            continue;
        }

        // 4. ピーク合わせ TGraph の構築
        vector<Double_t> x(_DT5751Length);
        vector<Double_t> y(_DT5751Length);
 
        for (Int_t k = 0; k < _DT5751Length; ++k) {
            x[k] = (Double_t)(k - k_peak); // 横軸: ピーク位置を 0 ns に揃える
            if (norm_flag) {
                y[k] = wave[k] / max_val; // 縦軸: ピーク高を 1.0 に規格化
            } else {
                y[k] = wave[k]; // 縦軸: 生のパルス高 (ADC)
            }
        }

        TGraph* g = new TGraph(_DT5751Length, &x[0], &y[0]);
        Int_t col = colors[count % colors.size()];
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
    
    // 積分ゲート範囲を示す縦線 (TLine) の描画
    Double_t ymin = norm_flag ? -0.1 : mg->GetYaxis()->GetXmin();
    Double_t ymax = norm_flag ? 1.1 : mg->GetYaxis()->GetXmax();

    auto draw_gate_line = [](Double_t x_pos, Double_t y_min, Double_t y_max, Int_t color) {
        TLine* line = new TLine(x_pos, y_min, x_pos, y_max);
        line->SetLineColor(color);
        line->SetLineWidth(2);
        line->SetLineStyle(2); // 破線
        line->Draw("same");
    };

    // 積分開始 (Peak - pre_ns) -> 赤色破線
    draw_gate_line(-(Double_t)n_pre_peak, ymin, ymax, kRed+1);
    // Short積分終了 (Peak + short_ns) -> 青色破線
    draw_gate_line((Double_t)n_post_peak_short, ymin, ymax, kBlue+1);
    // Long積分終了 (Peak + long_ns) -> 緑色破線
    draw_gate_line((Double_t)n_post_peak_long, ymin, ymax, kGreen+2);

    c->Update();
    
    // 入力ファイル名からベース名を取得 (例: data/Cf252_wave_01_slown.root -> Cf252_wave_01_slown)
    string base_name = input_path;
    size_t last_slash = base_name.find_last_of("/\\");
    if (last_slash != string::npos) {
        base_name = base_name.substr(last_slash + 1);
    }
    size_t last_dot = base_name.find_last_of(".");
    if (last_dot != string::npos) {
        base_name = base_name.substr(0, last_dot);
    }

    // 画像およびPDFとして自動保存
    string out_pdf = "pdf/" + base_name + "_overlay_ch" + to_string(target_ch) + (norm_flag ? "_norm" : "") + ".pdf";
    gSystem->mkdir("pdf", true);
    c->Print(out_pdf.c_str());
    cout << "Overlay plot saved to: " << out_pdf << endl;

    // ウインドウ表示を維持
    app.Run();

    file->Close();
    delete file;
    return 0;
}
