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
#include "progress_bar.h"

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

    Double_t delta_T_us_val = -1.0;
    bool has_delta_T = false;
    if (wave_tree->GetBranch("delta_T_us")) {
        wave_tree->SetBranchAddress("delta_T_us", &delta_T_us_val);
        has_delta_T = true;
    }

    bool is_slown = (filepath.find("slown") != string::npos);

    Int_t n_ranges = (Int_t)ceil(max_height / height_step);
    collector.bucket_waves.resize(n_ranges);

    Long64_t n_entries = wave_tree->GetEntries();
    const Double_t threshold = 5.0;
    const Int_t n_baseline_length = 200;

    cout << "\nScanning file: " << filepath << "..." << endl;
    for (Long64_t i = 0; i < n_entries; ++i) {
        if (i % 5000 == 0 || i == n_entries - 1) {
            displayProgressBar(i + 1, n_entries);
        }

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

        // delta_T_us のカット (slownでは10us未満の即発ガンマ背景をカット)
        if (has_delta_T) {
            if (is_slown) {
                if (delta_T_us_val < 10.0) continue;
            } else {
                if (delta_T_us_val < 0.0) continue;
            }
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

    string file1_path = "root/Co60_wave_merge_gamma.root";
    string file2_path = "root/Cf252_wave_merge_slown.root";
    Int_t target_ch = 1;
    Double_t height_step = 50.0;
    Double_t max_height = 1000.0;
    Int_t num_to_average = 1000000; // デフォルトで実質上限なし (使える波形は全部使う)

    Int_t opt_start = 1;
    if (argc >= 3 && argv[1][0] != '-' && argv[2][0] != '-') {
        file1_path = argv[1];
        file2_path = argv[2];
        opt_start = 3;
    } else if (argc >= 2 && argv[1][0] != '-') {
        file1_path = argv[1];
        opt_start = 2;
    }

    for (Int_t i = opt_start; i < argc; ++i) {
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

    // ファイルラベルの自動抽出と明示的な差し替え
    auto get_display_label = [](const string& path) -> string {
        size_t slash = path.find_last_of("/\\");
        string base = (slash == string::npos) ? path : path.substr(slash + 1);
        size_t dot = base.find_last_of(".");
        if (dot != string::npos) {
            base = base.substr(0, dot);
        }
        
        if (base.find("Co60") != string::npos || base.find("Co_60") != string::npos) {
            return "Co-60 gamma";
        } else if (base.find("Cf252") != string::npos) {
            if (base.find("slown") != string::npos) {
                return "Cf-252 slown";
            } else if (base.find("fastn") != string::npos) {
                return "Cf-252 fastn";
            } else if (base.find("gamma") != string::npos) {
                return "Cf-252 gamma";
            }
            return "Cf-252";
        }
        return base;
    };

    WaveformCollector col1;
    col1.file_label = get_display_label(file1_path);
    col1.color = kBlue; // ファイル1 (Gamma) -> 青色

    WaveformCollector col2;
    col2.file_label = get_display_label(file2_path);
    col2.color = kRed;  // ファイル2 (Slow Neutron) -> 赤色

    // 波形収集の実行
    if (!collect_waveforms(file1_path, target_ch, height_step, max_height, num_to_average, col1)) return 1;
    if (!collect_waveforms(file2_path, target_ch, height_step, max_height, num_to_average, col2)) return 1;

    Int_t n_ranges = (Int_t)ceil(max_height / height_step);

    // 出力PDFの保存先パス自動生成 (file1_pathのフォルダ構造を維持して root を pdf に置換)
    auto get_basename = [](const string& path) {
        size_t slash = path.find_last_of("/\\");
        string base = (slash == string::npos) ? path : path.substr(slash + 1);
        size_t dot = base.find_last_of(".");
        return (dot == string::npos) ? base : base.substr(0, dot);
    };

    string out_pdf;
    string base_name = file1_path;
    size_t last_dot = base_name.find_last_of(".");
    if (last_dot != string::npos) {
        base_name = base_name.substr(0, last_dot);
    }
    
    string suffix = "_vs_" + get_basename(file2_path) + ".pdf";
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
        out_pdf = "pdf/compare_average_" + base_name + suffix;
    }

    // 出力先フォルダの作成
    size_t last_slash_pdf = out_pdf.find_last_of("/\\");
    if (last_slash_pdf != string::npos) {
        string pdf_dir = out_pdf.substr(0, last_slash_pdf);
        gSystem->mkdir(pdf_dir.c_str(), true);
    } else {
        gSystem->mkdir("pdf", true);
    }

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
        // 凡例を左上のデッドスペース（立ち上がり前）に配置し、枠内に綺麗に収める
        TLegend* leg = new TLegend(0.15, 0.72, 0.48, 0.86);
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

        // 差分波形 (col2 - col1 = Cf - Co) の計算
        TGraphErrors* ge_diff = nullptr;
        if (ge1 && ge2) {
            Double_t* x1 = ge1->GetX();
            Double_t* y1 = ge1->GetY();
            Double_t* ey1 = ge1->GetEY();
            Double_t* y2 = ge2->GetY();
            Double_t* ey2 = ge2->GetEY();

            vector<Double_t> x_diff(_DT5751Length);
            vector<Double_t> y_diff(_DT5751Length);
            vector<Double_t> ex_diff(_DT5751Length, 0.0);
            vector<Double_t> ey_diff(_DT5751Length);

            Double_t diff_max = 0.0;
            for (Int_t idx = 0; idx < _DT5751Length; ++idx) {
                x_diff[idx] = x1[idx];
                y_diff[idx] = y2[idx] - y1[idx]; // Cf - Co の差分
                ey_diff[idx] = sqrt(pow(ey1[idx], 2) + pow(ey2[idx], 2)); // 誤差伝播
                if (y_diff[idx] > diff_max) {
                    diff_max = y_diff[idx];
                }
            }

            ge_diff = new TGraphErrors(_DT5751Length, &x_diff[0], &y_diff[0], &ex_diff[0], &ey_diff[0]);
            ge_diff->SetLineColor(kBlack); // 差分は黒色太線
            ge_diff->SetLineWidth(3);
            ge_diff->SetFillColorAlpha(kBlack, 0.15); // 差分のエラーバンド
            ge_diff->SetFillStyle(1001);

        }

        // ------------------------------------
        // ページ1: 重ね描き比較プロット
        // ------------------------------------
        mg->Draw("A");
        mg->GetXaxis()->SetRangeUser(-20.0, 100.0);
        mg->GetYaxis()->SetRangeUser(-20.0, range_max * 1.2); // Y軸のダイナミックレンジ調整

        leg->Draw("same");

        c->Update();
        c->Print(out_pdf.c_str());
        total_pages++;

        // ------------------------------------
        // ページ2: 差分単体プロット (差分データがある場合のみ)
        // ------------------------------------
        if (ge_diff) {
            c->Clear();
            c->SetGrid();

            TMultiGraph* mg_diff = new TMultiGraph();
            string diff_title;
            if (k == n_ranges - 1) {
                diff_title = Form("Difference (%s - %s) (Height: > %.0f ADC) [CH%d];Time relative to peak [ns];Difference [ADC]", 
                                  col2.file_label.c_str(), col1.file_label.c_str(), range_min, target_ch);
            } else {
                diff_title = Form("Difference (%s - %s) (Height: %.0f - %.0f ADC) [CH%d];Time relative to peak [ns];Difference [ADC]", 
                                  col2.file_label.c_str(), col1.file_label.c_str(), range_min, range_max, target_ch);
            }
            mg_diff->SetTitle(diff_title.c_str());
            mg_diff->Add(ge_diff, "3L");

            mg_diff->Draw("A");
            mg_diff->GetXaxis()->SetRangeUser(-20.0, 100.0);
            
            // 差分のピーク値に応じたダイナミックオートスケール (つぶれを完全に解消)
            Double_t y_max_diff = max(10.0, diff_max * 1.25);
            Double_t y_min_diff = -max(5.0, diff_max * 0.15);
            mg_diff->GetYaxis()->SetRangeUser(y_min_diff, y_max_diff);

            TLegend* leg_diff = new TLegend(0.15, 0.78, 0.48, 0.86);
            leg_diff->SetFillStyle(0);
            leg_diff->SetBorderSize(0);
            leg_diff->SetTextSize(0.032);
            leg_diff->AddEntry(ge_diff, "Cf - Co (Diff #pm #sigma)", "FL");
            leg_diff->Draw("same");

            c->Update();
            c->Print(out_pdf.c_str());
            total_pages++;

            delete mg_diff;
            delete leg_diff;
        }

        if (ge1) delete ge1;
        if (ge2) delete ge2;
        if (ge_diff) delete ge_diff;
    }

    c->Print(Form("%s)", out_pdf.c_str())); // PDFマルチページ完了

    cout << "Finished comparison. Generated " << total_pages << " pages." << endl;
    cout << "Output saved to: " << out_pdf << endl;

    delete c;
    return 0;
}
