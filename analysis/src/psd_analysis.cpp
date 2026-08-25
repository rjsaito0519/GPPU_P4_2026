#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

using namespace std;

static const int _DT5751Length = 1029;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_waveform_root_file> [output_root_file]" << endl;
        return 1;
    }

    string input_path = argv[1];
    string output_path = "";
    if (argc > 2) {
        output_path = argv[2];
    } else {
        output_path = input_path;
        size_t last_dot = output_path.find_last_of(".");
        if (last_dot != string::npos) {
            output_path = output_path.substr(0, last_dot);
        }
        output_path += "_psd.root";
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

    // パラメータ定義 (PSD.pdfおよびwave2tq.cc準拠)
    const double threshold = 10.0; // ADC value (threshold to identify real signal)
    const int n_baseline_length = 200;
    const double impedance = 50.0; // ohm
    const double dt = 1.0; // ns

    // ゲート定義（ピーク波高の割合で決定）
    const int n_pre_peak = 10;          // ピーク手前の積分開始オフセット
    const double fraction_short = 0.50; // Q_shortの終了点: ピーク高の50%以下まで減衰した位置
    const double fraction_long = 0.10;  // Q_longの終了点: ピーク高の10%以下まで減衰した位置

    Int_t event;
    ULong64_t time_stamp;
    Int_t channel;
    UShort_t wave_raw[_DT5751Length];

    wave_tree->SetBranchAddress("event", &event);
    wave_tree->SetBranchAddress("channel", &channel);
    wave_tree->SetBranchAddress("time_stamp", &time_stamp);
    wave_tree->SetBranchAddress("wave_raw", wave_raw);

    // 出力ROOTファイルのオープン
    TFile* outFile = TFile::Open(output_path.c_str(), "RECREATE");
    TTree* psd_tree = new TTree("psd_tree", "PSD Analysis Tree (CH1 only)");

    Double_t Q_long;
    Double_t Q_short;
    Double_t PSD;
    Double_t baseline;
    Double_t peak_time;

    psd_tree->Branch("event", &event, "event/I");
    psd_tree->Branch("time_stamp", &time_stamp, "time_stamp/l");
    psd_tree->Branch("Q_long", &Q_long, "Q_long/D");
    psd_tree->Branch("Q_short", &Q_short, "Q_short/D");
    psd_tree->Branch("PSD", &PSD, "PSD/D");
    psd_tree->Branch("baseline", &baseline, "baseline/D");
    psd_tree->Branch("peak_time", &peak_time, "peak_time/D");

    Long64_t n_entries = wave_tree->GetEntries();
    cout << "Analyzing waveforms (CH1 only) and calculating PSD..." << endl;
    cout << " - Short gate: " << gate_short << " ns" << endl;
    cout << " - Long gate: " << gate_long << " ns" << endl;
    cout << " - Target Output: " << output_path << endl;

    Long64_t analyzed_count = 0;

    for (Long64_t i = 0; i < n_entries; ++i) {
        wave_tree->GetEntry(i);

        // 制限条件: CH1のみ解析 (CH0や他のCHはスキップ)
        if (channel != 1) {
            continue;
        }

        // 1. ベースライン計算 (wave2tq.ccのロジック)
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
        baseline = (n_wave > 0) ? (sum_wave / double(n_wave)) : baseline_rough;

        // 2. ベースライン差分（極性反転）波形の作成
        vector<double> wave(_DT5751Length);
        for (int k = 0; k < _DT5751Length; ++k) {
            wave[k] = baseline - (double)wave_raw[k];
        }

        // 3. ピーク位置の探索 (400ns近辺 = 300 ~ 500 ns の範囲)
        double max_val = -99999.0;
        int k_peak = -1;
        for (int k = 300; k < 500; ++k) {
            if (wave[k] > max_val) {
                max_val = wave[k];
                k_peak = k;
            }
        }

        // 有意なパルスピークがない（ノイズ閾値以下）場合は解析対象外としてスキップ
        if (k_peak == -1 || max_val < threshold) {
            continue;
        }

        peak_time = (Double_t)k_peak;

        // 4. 積算電荷（Q_short, Q_long）の計算
        // 積分開始インデックス: ピーク手前 10 ns
        int k_start = k_peak - n_pre_peak;
        if (k_start < 0) k_start = 0;

        // Q_shortの積分終了点の決定 (ピーク以降で50%以下になる点)
        int k_short_end = k_peak;
        while (k_short_end < _DT5751Length && wave[k_short_end] > fraction_short * max_val) {
            k_short_end++;
        }

        // Q_longの積分終了点の決定 (ピーク以降で10%以下になる点)
        int k_long_end = k_peak;
        while (k_long_end < _DT5751Length && wave[k_long_end] > fraction_long * max_val) {
            k_long_end++;
        }

        Q_short = 0.0;
        Q_long = 0.0;

        // Q_shortの積分
        for (int k = k_start; k < k_short_end; ++k) {
            Q_short += (wave[k] / impedance * dt);
        }

        // Q_longの積分
        for (int k = k_start; k < k_long_end; ++k) {
            Q_long += (wave[k] / impedance * dt);
        }

        // 5. PSD パラメータ計算: (Q_long - Q_short) / Q_long
        if (Q_long > 0.0) {
            PSD = (Q_long - Q_short) / Q_long;
        } else {
            PSD = 0.0;
        }

        psd_tree->Fill();
        analyzed_count++;
    }

    psd_tree->Write();
    outFile->Close();
    file->Close();

    delete outFile;
    delete file;

    cout << "Finished PSD analysis. Analyzed " << analyzed_count << " CH1 events." << endl;
    cout << "Results saved to: " << output_path << endl;

    return 0;
}
