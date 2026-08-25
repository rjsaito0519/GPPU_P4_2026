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
        cerr << "Usage: " << argv[0] << " <input_waveform_root_file> [output_root_file] [--pre pre_ns] [--short short_ns] [--long long_ns] [--smooth 0|1]" << endl;
        return 1;
    }

    string input_path = argv[1];
    string output_path = "";
    Int_t n_pre_peak = 10;          // ピーク手前の積分開始オフセット [ns]
    Int_t n_post_peak_short = 10;   // ピーク後のQ_short積分幅 [ns]
    Int_t n_post_peak_long = 30;    // ピーク後のQ_long積分幅 [ns]
    Int_t apply_smoothing = 1;      // デジタル平滑化の有効化 (1: On, 0: Off)

    // オプション引数の賢いパース
    for (Int_t i = 2; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--pre" && i + 1 < argc) {
            n_pre_peak = stoi(argv[++i]);
        } else if (arg == "--short" && i + 1 < argc) {
            n_post_peak_short = stoi(argv[++i]);
        } else if (arg == "--long" && i + 1 < argc) {
            n_post_peak_long = stoi(argv[++i]);
        } else if (arg == "--smooth" && i + 1 < argc) {
            apply_smoothing = stoi(argv[++i]);
        } else {
            if (arg.find("--") == string::npos) {
                output_path = arg;
            }
        }
    }

    if (output_path.empty()) {
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

    TTree* wave_tree = (TTree*)file->Get("tree");
    if (!wave_tree) {
        cerr << "ERROR: cannot find TTree 'wave_tree' in input file" << endl;
        file->Close();
        return 1;
    }

    // パラメータ定義 (PSD.pdfおよびwave2tq.cc準拠)
    const Double_t threshold = 10.0; // ADC value (threshold to identify real signal)
    const Int_t n_baseline_length = 200;
    const Double_t impedance = 50.0; // ohm
    const Double_t dt = 1.0; // ns

    // wave2tq互換固定ゲートQ1用のパラメータ
    const Int_t n_softtrigger_enable_length = 800; // software trigger enable length
    const Int_t n_length_pre_softtrigger    =  10; // charge window length (pre)
    const Int_t n_length_post_softtrigger   = 100; // charge window length (post)

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
    TTree* psd_tree = new TTree("tree", "PSD Analysis Tree (CH1 only)");

    Double_t Q_long;
    Double_t Q_short;
    Double_t PSD;
    Double_t Q1;
    Double_t baseline;
    Double_t peak_time;
    Double_t t_short; // ピークからQ_short判定点までの時間差 (ns)
    Double_t t_long;  // ピークからQ_long判定点までの時間差 (ns)

    psd_tree->Branch("event", &event, "event/I");
    psd_tree->Branch("time_stamp", &time_stamp, "time_stamp/l");
    psd_tree->Branch("Q_long", &Q_long, "Q_long/D");
    psd_tree->Branch("Q_short", &Q_short, "Q_short/D");
    psd_tree->Branch("Q1", &Q1, "Q1/D");
    psd_tree->Branch("t_short", &t_short, "t_short/D");
    psd_tree->Branch("t_long", &t_long, "t_long/D");
    psd_tree->Branch("PSD", &PSD, "PSD/D");
    psd_tree->Branch("baseline", &baseline, "baseline/D");
    psd_tree->Branch("peak_time", &peak_time, "peak_time/D");

    Long64_t n_entries = wave_tree->GetEntries();
    cout << "Analyzing waveforms (CH1 only) and calculating PSD..." << endl;
    cout << " - Short gate: [Peak - " << n_pre_peak << " ns] to [Peak + " << n_post_peak_short << " ns]" << endl;
    cout << " - Long gate: [Peak - " << n_pre_peak << " ns] to [Peak + " << n_post_peak_long << " ns]" << endl;
    cout << " - Smoothing (Low-pass filter): " << (apply_smoothing ? "ON (1:2:1 Weighted Average)" : "OFF") << endl;
    cout << " - Target Output: " << output_path << endl;

    Long64_t analyzed_count = 0;

    for (Long64_t i = 0; i < n_entries; ++i) {
        wave_tree->GetEntry(i);

        // 制限条件: CH1のみ解析 (CH0や他のCHはスキップ)
        if (channel != 1) {
            continue;
        }

        // 1. ベースライン計算 (wave2tq.ccのロジック)
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
        baseline = (n_wave > 0) ? (sum_wave / Double_t(n_wave)) : baseline_rough;

        // 2. ベースライン差分（極性反転）波形の作成
        vector<Double_t> wave(_DT5751Length);
        for (Int_t k = 0; k < _DT5751Length; ++k) {
            wave[k] = baseline - (Double_t)wave_raw[k];
        }

        // 3点移動平均による平滑化（ローパスフィルタ）
        if (apply_smoothing) {
            vector<Double_t> temp_wave = wave;
            for (Int_t k = 1; k < _DT5751Length - 1; ++k) {
                wave[k] = (temp_wave[k - 1] + 2.0 * temp_wave[k] + temp_wave[k + 1]) / 4.0;
            }
        }

        // 3. ピーク位置の探索 (400ns近辺 = 300 ~ 500 ns の範囲)
        Double_t max_val = -99999.0;
        Int_t k_peak = -1;
        for (Int_t k = 300; k < 500; ++k) {
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
        Int_t k_start = k_peak - n_pre_peak;
        if (k_start < 0) k_start = 0;

        // Q_shortの積分終了点 (ピーク後 n_post_peak_short ns)
        Int_t k_short_end = min(_DT5751Length, k_peak + n_post_peak_short);
 
        // Q_longの積分終了点 (ピーク後 n_post_peak_long ns)
        Int_t k_long_end = min(_DT5751Length, k_peak + n_post_peak_long);
 
        // ピークトップから積分終了判定ポイントまでの時間差 (ns) の計算
        t_short = (Double_t)(k_short_end - k_peak) * dt;
        t_long = (Double_t)(k_long_end - k_peak) * dt;
 
        Q_short = 0.0;
        Q_long = 0.0;
 
        // Q_shortの積分
        for (Int_t k = k_start; k < k_short_end; ++k) {
            Q_short += (wave[k] / impedance * dt);
        }
 
        // Q_longの積分
        for (Int_t k = k_start; k < k_long_end; ++k) {
            Q_long += (wave[k] / impedance * dt);
        }

        // 4.1 wave2tq互換の固定窓電荷(Q1)の計算
        Q1 = 0.0;
        Int_t k_threshold = -1;
        for (Int_t k = 0; k < _DT5751Length; ++k) {
            if (k >= n_softtrigger_enable_length) break;
            if (wave[k] >= threshold) {
                k_threshold = k;
                break;
            }
        }
        if (k_threshold != -1) {
            Int_t k_window_start = k_threshold - n_length_pre_softtrigger;
            Int_t k_window_end = k_threshold + n_length_post_softtrigger;
            if (k_window_start < 0) k_window_start = 0;
            if (k_window_end > _DT5751Length) k_window_end = _DT5751Length;
            
            for (Int_t k = k_window_start; k < k_window_end; ++k) {
                Q1 += (wave[k] / impedance * dt);
            }
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
