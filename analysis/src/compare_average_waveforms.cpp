#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TLegend.h>

using namespace std;

static const Int_t _DT5751Length = 1029;

// 各ファイルから集められた波形を格納する構造体
struct WaveformCollector {
    string file_label;
    Int_t color;
    vector<vector<vector<Double_t>>> bucket_waves; // [range_idx][event_idx][sample_idx]
};

// 指定ファイルから波高レンジごとに波形を収集する関数
bool collect_waveforms(const string& filepath, Int_t target_ch, Double_t height_step, Double_t max_height, Int_t num_to_average, WaveformCollector& collector) {
    TFile* file = TFile::Open(filepath.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "ERROR: cannot open file -> " << filepath << endl;
        return false;
    }

    TTree* wave_tree = (TTree*)file->Get("tree");
    if (!wave_tree) {
        cerr << "ERROR: cannot find TTree 'tree' in file -> " << filepath << endl;
        file->Close();
        return false;
    }

    Int_t event;
    Int_t channel;
    ULong64_t time_stamp;
    UShort_t wave_raw[_DT5751Length];

    wave_tree->SetBranchAddress("event", &event);
    wave_tree->SetBranchAddress("channel", &channel);
    wave_tree->SetBranchAddress("time_stamp", &time_stamp);
    wave_tree->SetBranchAddress("wave_raw", wave_raw);

    Int_t n_ranges = (Int_t)ceil(max_height / height_step);
    collector.bucket_waves.resize(n_ranges);

    Long64_t n_entries = wave_tree->GetEntries();
    const Double_t threshold = 5.0;
    const Int_t n_baseline_length = 200;

    for (Long64_t i = 0; i < n_entries; ++i) {
        // すべてのバケツが満杯か確認
        bool all_full = true;
        for (Int_t k = 0; k < n_ranges; ++k) {
            if ((Int_t)collector.bucket_waves[k].size() < num_to_average) {
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
 
        sum_wave = 0.0; n_wave = 0;
        for (Int_t k = 0; k < n_baseline_length; ++k) {
            if (abs((Double_t)wave_raw[k] - baseline_rough) < threshold) {
                sum_wave += (Double_t)wave_raw[k];
                n_wave++;
            }
        }
        Double_t baseline = (n_wave > 0) ? (sum_wave / Double_t(n_wave)) : baseline_rough;
 
        // 2. 反転差分波形
        vector<Double_t> wave(_DT5751Length);
        for (Int_t k = 0; k < _DT5751Length; ++k) {
            wave[k] = baseline - (Double_t)wave_raw[k];
        }
 
        // 3. ピークサーチ
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

        // 対応バケットインデックスの決定
        Int_t range_idx = (Int_t)(max_val / height_step);
        if (range_idx >= n_ranges) {
            range_idx = n_ranges - 1; // 1000以上の飽和パルスは最後のバケツへ
        }

        if (range_idx >= 0 && (Int_t)collector.bucket_waves[range_idx].size() < num_to_average) {
            // アライメントを考慮した波形データの作成
            vector<Double_t> aligned_wave(_DT5751Length, 0.0);
            for (Int_t k = 0; k < _DT5751Length; ++k) {
                // ピークを k = 200 にシフトして時間アライメントを取る
                Int_t orig_idx = k - 200 + k_peak;
                if (orig_idx >= 0 && orig_idx < _DT5751Length) {
                    aligned_wave[k] = wave[orig_idx];
                }
            }
            collector.bucket_waves[range_idx].push_back(aligned_wave);
        }
    }

    file->Close();
    delete file;
    return true;
}

// 平均波形および標準偏差エラーのグラフを生成する関数
TGraphErrors* calculate_average_graph(const vector<vector<Double_t>>& waves, Int_t color) {
    if (waves.empty()) return nullptr;

    Int_t n_events = waves.size();
    vector<Double_t> x(_DT5751Length);
    vector<Double_t> y(_DT5751Length, 0.0);
    vector<Double_t> ex(_DT5751Length, 0.0);
    vector<Double_t> ey(_DT5751Length, 0.0);

    for (Int_t k = 0; k < _DT5751Length; ++k) {
        x[k] = (Double_t)(k - 200); // ピーク位置 (k=200) を 0 ns に設定

        // 平均値 (Mean) の計算
        Double_t sum = 0.0;
        for (Int_t i = 0; i < n_events; ++i) {
            sum += waves[i][k];
        }
        y[k] = sum / (Double_t)n_events;

        // 標準偏差 (Standard Deviation) の計算 (統計誤差表現)
        if (n_events > 1) {
            Double_t var_sum = 0.0;
            for (Int_t i = 0; i < n_events; ++i) {
                var_sum += pow(waves[i][k] - y[k], 2);
            }
            ey[k] = sqrt(var_sum / (Double_t)(n_events - 1));
        } else {
            ey[k] = 0.0;
        }
    }

    TGraphErrors* ge = new TGraphErrors(_DT5751Length, &x[0], &y[0], &ex[0], &ey[0]);
    ge->SetLineColor(color);
    ge->SetLineWidth(3);
    ge->SetFillColorAlpha(color, 0.25); // エラーバンド用透過色設定
    ge->SetFillStyle(1001);
    return ge;
}

int main(int argc, char* argv[]) {
    gROOT->SetBatch(kTRUE);

    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <file1_root> <file2_root> [--ch channel] [--step height_step] [--max max_height] [--n num_to_average]" << endl;
        return 1;
    }

    string file1_path = argv[1];
    string file2_path = argv[2];
    Int_t target_ch = 1;
    Double_t height_step = 50.0;
    Double_t max_height = 1000.0;
    Int_t num_to_average = 100; // 平均値算出用のパルス上限数 (多いほどノイズが消えます)

    for (Int_t i = 3; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--ch" && i + 1 < argc) {
            target_ch = stoi(argv[++i]);
        } else if (arg == "--step" && i + 1 < argc) {
            height_step = stod(argv[++i]);
        } else if (arg == "--max" && i + 1 < argc) {
            max_height = stod(argv[++i]);
        } else if (arg == "--n" && i + 1 < argc) {
            num_to_average = stoi(argv[++i]);
        }
    }

    // ファイルラベルの自動抽出
    auto get_basename = [](const string& path) {
        size_t slash = path.find_last_of("/\\");
        string base = (slash == string::npos) ? path : path.substr(slash + 1);
        size_t dot = base.find_last_of(".");
        return (dot == string::npos) ? base : base.substr(0, dot);
    };

    WaveformCollector col1;
    col1.file_label = get_basename(file1_path);
    col1.color = kBlue; // ファイル1 (Gamma) -> 青色

    WaveformCollector col2;
    col2.file_label = get_basename(file2_path);
    col2.color = kRed;  // ファイル2 (Slow Neutron) -> 赤色

    // 波形収集の実行
    if (!collect_waveforms(file1_path, target_ch, height_step, max_height, num_to_average, col1)) return 1;
    if (!collect_waveforms(file2_path, target_ch, height_step, max_height, num_to_average, col2)) return 1;

    Int_t n_ranges = (Int_t)ceil(max_height / height_step);

    // 出力PDFパス決定 (pdf/others/ に自動整理)
    string out_pdf = "pdf/others/compare_average_" + col1.file_label + "_vs_" + col2.file_label + ".pdf";
    gSystem->mkdir("pdf/others", true);

    TCanvas* c = new TCanvas("c_compare", "Average Waveform Comparison", 900, 700);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetGrid();

    cout << "Comparing average waveforms with error bands..." << endl;
    c->Print(Form("%s(", out_pdf.c_str())); // PDFマルチページ開始

    Int_t total_pages = 0;
    for (Int_t k = 0; k < n_ranges; ++k) {
        if (col1.bucket_waves[k].empty() && col2.bucket_waves[k].empty()) {
            continue; // 両ファイルともデータがないレンジはスキップ
        }

        c->Clear();
        c->SetGrid();

        TMultiGraph* mg = new TMultiGraph();
        TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.032);

        Double_t range_min = k * height_step;
        Double_t range_max = (k == n_ranges - 1) ? max_height : ((k + 1) * height_step);

        string title;
        if (k == n_ranges - 1) {
            title = Form("Average Pulse (Height: > %.0f ADC) [CH%d];Time relative to peak [ns];Pulse Height [ADC]", range_min, target_ch);
        } else {
            title = Form("Average Pulse (Height: %.0f - %.0f ADC) [CH%d];Time relative to peak [ns];Pulse Height [ADC]", range_min, range_max, target_ch);
        }
        mg->SetTitle(title.c_str());

        // ファイル1の平均波形とエラーバンドの計算＆追加
        TGraphErrors* ge1 = calculate_average_graph(col1.bucket_waves[k], col1.color);
        if (ge1) {
            mg->Add(ge1, "3L"); // 3: エラーバンド, L: 平均線
            leg->AddEntry(ge1, Form("%s (N=%d)", col1.file_label.c_str(), (int)col1.bucket_waves[k].size()), "FL");
        }

        // ファイル2の平均波形とエラーバンドの計算＆追加
        TGraphErrors* ge2 = calculate_average_graph(col2.bucket_waves[k], col2.color);
        if (ge2) {
            mg->Add(ge2, "3L");
            leg->AddEntry(ge2, Form("%s (N=%d)", col2.file_label.c_str(), (int)col2.bucket_waves[k].size()), "FL");
        }

        mg->Draw("A");
        mg->GetXaxis()->SetRangeUser(-20.0, 150.0);
        mg->GetYaxis()->SetRangeUser(-20.0, range_max * 1.2); // Y軸のダイナミックレンジ調整

        leg->Draw("same");

        c->Update();
        c->Print(out_pdf.c_str());
        total_pages++;

        if (ge1) delete ge1;
        if (ge2) delete ge2;
    }

    c->Print(Form("%s)", out_pdf.c_str())); // PDFマルチページ完了

    cout << "Finished comparison. Generated " << total_pages << " pages." << endl;
    cout << "Output saved to: " << out_pdf << endl;

    delete c;
    return 0;
}
