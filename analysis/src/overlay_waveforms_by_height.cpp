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
#include <TAxis.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>

using namespace std;

static const Int_t _DT5751Length = 1029;
 
int main(int argc, char* argv[]) {
    // X11ウインドウを開かないバッチモードを有効化
    gROOT->SetBatch(kTRUE);

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_root_file> [--n num_events] [--ch channel] [--norm 0|1] [--step height_step] [--max max_height]" << endl;
        return 1;
    }
 
    string input_path = argv[1];
    Int_t num_to_overlay = 20;       // 各レンジで重ね合わせるパルス数
    Int_t target_ch = 1;             // 対象チャンネル (液体シンチ)
    Int_t norm_flag = 0;             // デフォルトで波高比較のため規格化は OFF (0: Raw, 1: Normalized)
    Double_t height_step = 50.0;     // 波高の刻み幅 (デフォルト: 50 ADC)
    Double_t max_height = 1000.0;    // スキャンの最大波高制限 (デフォルト: 1000 ADC)
 
    for (Int_t i = 2; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--n" && i + 1 < argc) {
            num_to_overlay = stoi(argv[++i]);
        } else if (arg == "--ch" && i + 1 < argc) {
            target_ch = stoi(argv[++i]);
        } else if (arg == "--norm" && i + 1 < argc) {
            norm_flag = stoi(argv[++i]);
        } else if (arg == "--step" && i + 1 < argc) {
            height_step = stod(argv[++i]);
        } else if (arg == "--max" && i + 1 < argc) {
            max_height = stod(argv[++i]);
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

    Int_t event;
    Int_t channel;
    ULong64_t time_stamp;
    UShort_t wave_raw[_DT5751Length];

    wave_tree->SetBranchAddress("event", &event);
    wave_tree->SetBranchAddress("channel", &channel);
    wave_tree->SetBranchAddress("time_stamp", &time_stamp);
    wave_tree->SetBranchAddress("wave_raw", wave_raw);

    // バケット数 (レンジ数) の決定
    Int_t n_ranges = (Int_t)ceil(max_height / height_step);
    if (n_ranges <= 0) n_ranges = 1;

    // 波高レンジごとの TGraph 格納用バケット
    vector<vector<TGraph*>> bucket_graphs(n_ranges);

    Long64_t n_entries = wave_tree->GetEntries();
    const Double_t threshold = 5.0; // 有意パルスの判定閾値 (ADC)
    const Int_t n_baseline_length = 200;

    // カラーパレット (重ね描き用)
    vector<Int_t> colors = {kBlue, kRed, kGreen+2, kOrange+1, kMagenta, kCyan+2, kViolet, kPink-3, kTeal+9, kAzure+2};

    cout << "Scanning events and categorizing pulses by height (step: " << height_step << " ADC)..." << endl;

    for (Long64_t i = 0; i < n_entries; ++i) {
        // すべてのバケツが満杯に達したか確認
        bool all_full = true;
        for (Int_t k = 0; k < n_ranges; ++k) {
            if ((Int_t)bucket_graphs[k].size() < num_to_overlay) {
                all_full = false;
                break;
            }
        }
        if (all_full) break;

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

        // 対応する波高レンジバケットの判定
        Int_t range_idx = (Int_t)(max_val / height_step);
        if (range_idx >= n_ranges) {
            // max_heightを超える巨大パルスは最後のバケツに集約
            range_idx = n_ranges - 1;
        }

        if (range_idx >= 0 && (Int_t)bucket_graphs[range_idx].size() < num_to_overlay) {
            // TGraph の構築
            vector<Double_t> x(_DT5751Length);
            vector<Double_t> y(_DT5751Length);
     
            for (Int_t k = 0; k < _DT5751Length; ++k) {
                x[k] = (Double_t)(k - k_peak); // 横軸: ピーク位置を 0 ns に揃える
                if (norm_flag) {
                    y[k] = wave[k] / max_val; // 規格化
                } else {
                    y[k] = wave[k]; // 生の波高値
                }
            }

            TGraph* g = new TGraph(_DT5751Length, &x[0], &y[0]);
            Int_t color_idx = bucket_graphs[range_idx].size();
            g->SetLineColor(colors[color_idx % colors.size()]);
            g->SetLineWidth(1);

            bucket_graphs[range_idx].push_back(g);
        }
    }

    // PDF出力先パス自動生成 (data または root を pdf に自動置換して階層維持)
    string out_pdf;
    string base_name = input_path;
    size_t last_dot = base_name.find_last_of(".");
    if (last_dot != string::npos) {
        base_name = base_name.substr(0, last_dot);
    }
    
    string suffix = "_overlay_by_height_ch" + to_string(target_ch) + (norm_flag ? "_norm" : "") + ".pdf";

    size_t data_pos = base_name.find("data");
    size_t root_pos = base_name.find("root");
    if (data_pos != string::npos) {
        base_name.replace(data_pos, 4, "pdf");
        out_pdf = base_name + suffix;
    } else if (root_pos != string::npos) {
        base_name.replace(root_pos, 4, "pdf");
        out_pdf = base_name + suffix;
    } else {
        size_t last_slash = base_name.find_last_of("/\\");
        if (last_slash != string::npos) {
            base_name = base_name.substr(last_slash + 1);
        }
        out_pdf = "pdf/" + base_name + suffix;
    }

    // ディレクトリの作成
    size_t last_slash_pdf = out_pdf.find_last_of("/\\");
    if (last_slash_pdf != string::npos) {
        string pdf_dir = out_pdf.substr(0, last_slash_pdf);
        gSystem->mkdir(pdf_dir.c_str(), true);
    } else {
        gSystem->mkdir("pdf", true);
    }

    TCanvas* c = new TCanvas("c_height", "Waveform Overlay by Height", 900, 700);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetGrid();

    cout << "Generating PDF with overlaid pages..." << endl;

    // マルチページ PDF の書き出し開始
    c->Print(Form("%s(", out_pdf.c_str()));

    Int_t total_pages = 0;
    for (Int_t k = 0; k < n_ranges; ++k) {
        if (bucket_graphs[k].empty()) {
            continue; // 波形がないレンジはページを作成しない
        }

        c->Clear();
        c->SetGrid();

        TMultiGraph* mg = new TMultiGraph();
        
        Double_t range_min = k * height_step;
        Double_t range_max = (k == n_ranges - 1) ? (max_height * 5.0) : ((k + 1) * height_step);

        // 入力パスからベース名 (ファイル名のみ) を抽出
        string file_basename = input_path;
        size_t last_slash_in = file_basename.find_last_of("/\\");
        if (last_slash_in != string::npos) {
            file_basename = file_basename.substr(last_slash_in + 1);
        }

        string title;
        if (k == n_ranges - 1) {
            title = Form("%s: > %.0f ADC (CH%d)", file_basename.c_str(), range_min, target_ch);
        } else {
            title = Form("%s: %.0f - %.0f ADC (CH%d)", file_basename.c_str(), range_min, (k + 1) * height_step, target_ch);
        }
        if (norm_flag) title += " [Normalized]";
        title += ";Time relative to peak [ns];" + string(norm_flag ? "Normalized Amplitude" : "Pulse Height [ADC]");

        mg->SetTitle(title.c_str());

        for (auto g : bucket_graphs[k]) {
            mg->Add(g, "L");
        }

        mg->Draw("A");
        mg->GetXaxis()->SetRangeUser(-20.0, 150.0);

        if (norm_flag) {
            mg->GetYaxis()->SetRangeUser(-0.1, 1.1);
        } else {
            // Y軸の最大値制限を各ページの波高上限 (range_max) に追従させて見やすくする
            mg->GetYaxis()->SetRangeUser(-20.0, range_max * 1.2);
        }

        c->Update();
        c->Print(out_pdf.c_str());
        total_pages++;

        // 後片付け
        for (auto g : bucket_graphs[k]) {
            delete g;
        }
    }

    // マルチページ PDF の書き出し完了
    c->Print(Form("%s)", out_pdf.c_str()));

    cout << "Finished. Generated " << total_pages << " pages of overlaid plots." << endl;
    cout << "Output saved: " << out_pdf << endl;

    file->Close();
    delete file;
    delete c;

    return 0;
}
