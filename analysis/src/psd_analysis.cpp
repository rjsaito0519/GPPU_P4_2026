#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "progress_bar.h"

using namespace std;

static const int _DT5751Length = 1029;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_waveform_root_file> [output_root_file] [--pre pre_ns] [--short short_ns] [--long long_ns] [--smooth 0|1] [--quiet] [-n max_events]" << endl;
        return 1;
    }

    string input_path = argv[1];
    string output_path = "";
    Int_t n_pre_peak = 10;          // ピーク手前の積分開始オフセット [ns]
    Int_t n_post_peak_short = 10;   // ピーク後のQ_short積分幅 [ns]
    Int_t n_post_peak_long = 150;   // ピーク後のQ_long積分幅 [ns]
    Int_t apply_smoothing = 1;      // デジタル平滑化の有効化 (1: On, 0: Off)
    Bool_t quiet = false;           // プログレスバーの非表示フラグ
    Long64_t max_events = -1;       // 最大解析イベント数 (-1: 全件)

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
        } else if ((arg == "-n" || arg == "--n") && i + 1 < argc) {
            max_events = stoll(argv[++i]);
        } else if (arg == "--quiet") {
            quiet = true;
        } else {
            if (arg.find("--") == string::npos) {
                output_path = arg;
            }
        }
    }

    if (output_path.empty()) {
        string base_filename = input_path;
        size_t last_dot = base_filename.find_last_of(".");
        if (last_dot != string::npos) {
            base_filename = base_filename.substr(0, last_dot);
        }
        output_path = base_filename + "_psd.root";
    }

    // 出力先フォルダの作成
    size_t last_slash_out = output_path.find_last_of("/\\");
    if (last_slash_out != string::npos) {
        string out_dir = output_path.substr(0, last_slash_out);
        gSystem->mkdir(out_dir.c_str(), true);
    } else {
        gSystem->mkdir("root", true);
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

    // パラメータ定義 (PSD.pdfおよびwave2tq.cc準拠)
    const Double_t threshold = 10.0; // ADC value (threshold to identify real signal)
    const Int_t n_baseline_length = 200;
    const Double_t impedance = 50.0; // ohm
    const Double_t dt = 1.0; // ns

    // wave2tq互換固定ゲートQ用のパラメータ
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

    // 出力ROOTファイルのオープンと出力TTree定義
    TFile* outFile = TFile::Open(output_path.c_str(), "RECREATE");
    TTree* psd_tree = new TTree("tree", "TQ and PSD Composite Tree");

    // 出力変数 (グローバルバインド用)
    Int_t out_event = -1;
    ULong64_t out_time_stamp = 0;
    Double_t T0 = 0.0, Q0 = 0.0;
    Double_t T1, Q1;
    Double_t Q_long;
    Double_t Q_short;
    Double_t PSD;
    Double_t baseline;
    Double_t peak_time;
    Double_t t_short;
    Double_t t_long;

    psd_tree->Branch("event", &out_event, "event/I");
    psd_tree->Branch("time_stamp", &out_time_stamp, "time_stamp/l");
    psd_tree->Branch("T0", &T0, "T0/D");
    psd_tree->Branch("Q0", &Q0, "Q0/D");
    psd_tree->Branch("T1", &T1, "T1/D");
    psd_tree->Branch("Q1", &Q1, "Q1/D");
    psd_tree->Branch("Q_long", &Q_long, "Q_long/D");
    psd_tree->Branch("Q_short", &Q_short, "Q_short/D");
    psd_tree->Branch("t_short", &t_short, "t_short/D");
    psd_tree->Branch("t_long", &t_long, "t_long/D");
    psd_tree->Branch("PSD", &PSD, "PSD/D");
    psd_tree->Branch("baseline", &baseline, "baseline/D");
    psd_tree->Branch("peak_time", &peak_time, "peak_time/D");

    Long64_t n_entries = wave_tree->GetEntries();
    cout << "Analyzing waveforms (Stream processing)..." << endl;
    cout << " - Short gate: [Peak - " << n_pre_peak << " ns] to [Peak + " << n_post_peak_short << " ns]" << endl;
    cout << " - Long gate: [Peak - " << n_pre_peak << " ns] to [Peak + " << n_post_peak_long << " ns]" << endl;
    cout << " - Smoothing (Low-pass filter): " << (apply_smoothing ? "ON (1:2:1 Weighted Average)" : "OFF") << endl;
    if (max_events > 0) {
        cout << " - Max Events Limit: " << max_events << " (Total file entries: " << n_entries << ")" << endl;
    }
    cout << " - Target Output: " << output_path << endl;

    Long64_t analyzed_events = 0;
    Int_t current_event = -1;
    bool has_ch0 = false;
    bool has_ch1 = false;

    // 前のイベント結果を出力ツリーに詰め込んで状態をリセットするヘルパー関数
    auto fill_current_event = [&]() {
        if (current_event != -1 && (has_ch0 || has_ch1)) {
            out_event = current_event;
            psd_tree->Fill();
            analyzed_events++;
        }
        // 出力変数のクリア
        T0 = 0.0; Q0 = 0.0;
        T1 = 0.0; Q1 = 0.0;
        Q_long = 0.0; Q_short = 0.0; PSD = 0.0;
        t_short = 0.0; t_long = 0.0;
        baseline = 0.0; peak_time = 0.0;
        has_ch0 = false;
        has_ch1 = false;
    };

    for (Long64_t i = 0; i < n_entries; ++i) {
        // 最大イベント数に達した場合はループを抜ける
        if (max_events > 0 && analyzed_events >= max_events) {
            break;
        }

        if (!quiet && (i % 5000 == 0 || i == n_entries - 1)) {
            displayProgressBar(i + 1, n_entries);
        }

        wave_tree->GetEntry(i);

        // イベントIDが切り替わったら、これまでに集約したデータをFill
        if (event != current_event) {
            fill_current_event();
            current_event = event;
            out_time_stamp = time_stamp;
        }

        // --- A. CH0 (プラスチックトリガー等) の解析 ---
        if (channel == 0) {
            // ベースライン計算
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
            Double_t bl_ch0 = (n_wave > 0) ? (sum_wave / Double_t(n_wave)) : baseline_rough;

            // 反転ベースライン減算波形
            vector<Double_t> wave_ch0(_DT5751Length);
            for (Int_t k = 0; k < _DT5751Length; ++k) {
                wave_ch0[k] = bl_ch0 - (Double_t)wave_raw[k];
            }

            // しきい値を超えたタイミング探索
            Int_t k_th = -1;
            for (Int_t k = 0; k < n_softtrigger_enable_length; ++k) {
                if (wave_ch0[k] >= threshold) {
                    k_th = k;
                    break;
                }
            }

            if (k_th != -1) {
                // 固定ゲート積算 Q0 の計算
                Int_t k_w_start = max(0, k_th - n_length_pre_softtrigger);
                Int_t k_w_end = min(_DT5751Length, k_th + n_length_post_softtrigger);
                for (Int_t k = k_w_start; k < k_w_end; ++k) {
                    Q0 += (wave_ch0[k] / impedance * dt);
                }

                // 線形補間 T0 の計算
                if (k_th > 0) {
                    Double_t time_1 = (Double_t)(k_th - 1);
                    Double_t time_2 = (Double_t)k_th;
                    Double_t wave_1 = wave_ch0[k_th - 1];
                    Double_t wave_2 = wave_ch0[k_th];
                    Double_t slope = (wave_1 - wave_2) / (time_1 - time_2);
                    Double_t offset = (time_1 * wave_2 - time_2 * wave_1) / (time_1 - time_2);
                    if (slope != 0.0) {
                        T0 = (threshold - offset) / slope;
                    }
                }
            }
            has_ch0 = true;
        }

        // --- B. CH1 (液体シンチレータ) の解析 ---
        else if (channel == 1) {
            // ベースライン計算
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
            baseline = (n_wave > 0) ? (sum_wave / Double_t(n_wave)) : baseline_rough;

            // 反転ベースライン減算波形
            vector<Double_t> wave_ch1(_DT5751Length);
            for (Int_t k = 0; k < _DT5751Length; ++k) {
                wave_ch1[k] = baseline - (Double_t)wave_raw[k];
            }

            // しきい値を超えたタイミング探索 (T1)
            Int_t k_th = -1;
            for (Int_t k = 0; k < n_softtrigger_enable_length; ++k) {
                if (wave_ch1[k] >= threshold) {
                    k_th = k;
                    break;
                }
            }

            if (k_th != -1) {
                // 固定ゲート積算 Q1 の計算
                Int_t k_w_start = max(0, k_th - n_length_pre_softtrigger);
                Int_t k_w_end = min(_DT5751Length, k_th + n_length_post_softtrigger);
                for (Int_t k = k_w_start; k < k_w_end; ++k) {
                    Q1 += (wave_ch1[k] / impedance * dt);
                }

                // 線形補間 T1 の計算
                if (k_th > 0) {
                    Double_t time_1 = (Double_t)(k_th - 1);
                    Double_t time_2 = (Double_t)k_th;
                    Double_t wave_1 = wave_ch1[k_th - 1];
                    Double_t wave_2 = wave_ch1[k_th];
                    Double_t slope = (wave_1 - wave_2) / (time_1 - time_2);
                    Double_t offset = (time_1 * wave_2 - time_2 * wave_1) / (time_1 - time_2);
                    if (slope != 0.0) {
                        T1 = (threshold - offset) / slope;
                    }
                }
            }

            // --- 3. PSD 積分用波形処理（デジタル平滑化） ---
            vector<Double_t> wave_smooth = wave_ch1;
            if (apply_smoothing) {
                for (Int_t k = 1; k < _DT5751Length - 1; ++k) {
                    wave_smooth[k] = (wave_ch1[k - 1] + 2.0 * wave_ch1[k] + wave_ch1[k + 1]) / 4.0;
                }
            }

            // ピーク位置探索 (300 ~ 500 ns)
            Double_t max_val = -99999.0;
            Int_t k_peak = -1;
            for (Int_t k = 300; k < 500; ++k) {
                if (wave_smooth[k] > max_val) {
                    max_val = wave_smooth[k];
                    k_peak = k;
                }
            }

            // ピークが有意な場合のみ PSD 計算
            if (k_peak != -1 && max_val >= threshold) {
                peak_time = (Double_t)k_peak;

                Int_t k_start = max(0, k_peak - n_pre_peak);
                Int_t k_short_end = min(_DT5751Length, k_peak + n_post_peak_short);
                Int_t k_long_end = min(_DT5751Length, k_peak + n_post_peak_long);

                t_short = (Double_t)(k_short_end - k_peak) * dt;
                t_long = (Double_t)(k_long_end - k_peak) * dt;

                // Q_short/Q_long の平滑化電荷積算
                Double_t q_short_sum = 0.0;
                Double_t q_long_sum = 0.0;
                for (Int_t k = k_start; k < k_short_end; ++k) {
                    q_short_sum += (wave_smooth[k] / impedance * dt);
                }
                for (Int_t k = k_start; k < k_long_end; ++k) {
                    q_long_sum += (wave_smooth[k] / impedance * dt);
                }

                Q_short = q_short_sum;
                Q_long = q_long_sum;

                if (Q_long > 0.0) {
                    PSD = (Q_long - Q_short) / Q_long;
                }
            }
            has_ch1 = true;
        }
    }

    // 最後のイベント結果を掃き出し
    fill_current_event();

    psd_tree->Write();
    outFile->Close();
    file->Close();

    delete outFile;
    delete file;

    cout << "\nAnalysis complete. Output saved: " << output_path << endl;
    return 0;
}
