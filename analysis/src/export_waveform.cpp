#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

using namespace std;

static const Int_t _DT5751DataSize = 2090;
static const Int_t _DT5751Length = 1029;

typedef struct {
    unsigned int header[8];
    unsigned short waveform[_DT5751Length];
} DT5751WFdata;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_list_filename> [event_list_txt] [output_root_filename]" << endl;
        return 1;
    }

    string input_list_filename = argv[1];
    string event_list_txt = "";
    string target_mode = ""; // "n" や "gamma"
    string output_root_filename = "";

    // 引数の賢い判定
    if (argc == 3) {
        string arg2 = argv[2];
        if (arg2.size() > 4 && arg2.substr(arg2.size() - 4) == ".txt") {
            event_list_txt = arg2;
        } else if (arg2 == "n" || arg2 == "gamma" || arg2 == "g" || arg2 == "N" || arg2 == "GAMMA" ||
                   arg2 == "fastn" || arg2 == "fn" || arg2 == "FASTN" || arg2 == "FN" ||
                   arg2 == "slown" || arg2 == "sn" || arg2 == "SLOWN" || arg2 == "SN") {
            target_mode = arg2;
            if (target_mode == "N") target_mode = "n";
            if (target_mode == "GAMMA" || target_mode == "g") target_mode = "gamma";
            if (target_mode == "FASTN" || target_mode == "fn" || target_mode == "FN") target_mode = "fastn";
            if (target_mode == "SLOWN" || target_mode == "sn" || target_mode == "SN") target_mode = "slown";
        } else {
            output_root_filename = arg2;
        }
    } else if (argc > 3) {
        string arg2 = argv[2];
        if (arg2.size() > 4 && arg2.substr(arg2.size() - 4) == ".txt") {
            event_list_txt = arg2;
        } else {
            target_mode = arg2;
            if (target_mode == "N") target_mode = "n";
            if (target_mode == "GAMMA" || target_mode == "g") target_mode = "gamma";
            if (target_mode == "FASTN" || target_mode == "fn" || target_mode == "FN") target_mode = "fastn";
            if (target_mode == "SLOWN" || target_mode == "sn" || target_mode == "SN") target_mode = "slown";
        }
        output_root_filename = argv[3];
    }

    // サフィックスの自動抽出 (例: data/Cf252_tq_merge_n_events.txt -> suffix = "_n")
    string suffix = "";
    if (!target_mode.empty()) {
        suffix = "_" + target_mode;
    } else if (!event_list_txt.empty()) {
        size_t ev_pos = event_list_txt.find("_events.txt");
        if (ev_pos != string::npos) {
            size_t last_under = event_list_txt.find_last_of("_", ev_pos - 1);
            if (last_under != string::npos) {
                suffix = event_list_txt.substr(last_under, ev_pos - last_under); // "_n" や "_gamma" などを抽出
            }
        } else {
            size_t dot_pos = event_list_txt.find_last_of(".");
            size_t last_under = event_list_txt.find_last_of("_", dot_pos - 1);
            if (last_under != string::npos && last_under < dot_pos) {
                suffix = event_list_txt.substr(last_under, dot_pos - last_under);
            }
        }
    }

    // 出力先ROOTファイル名の自動生成 (サフィックスを付与し、root/ フォルダ直下に配置)
    if (output_root_filename.empty()) {
        string base_filename = input_list_filename;
        size_t last_slash_in = base_filename.find_last_of("/\\");
        if (last_slash_in != string::npos) {
            base_filename = base_filename.substr(last_slash_in + 1);
        }
        size_t last_dot = base_filename.find_last_of(".");
        if (last_dot != string::npos) {
            base_filename = base_filename.substr(0, last_dot);
        }
        output_root_filename = "root/" + base_filename + suffix + ".root";
    }

    // 出力先フォルダの作成
    gSystem->mkdir("root", true);

    string list_dir = "";
    size_t last_slash = input_list_filename.find_last_of("/\\");
    if (last_slash != string::npos) {
        list_dir = input_list_filename.substr(0, last_slash + 1); // "data/" など
    }

    ifstream flist(input_list_filename.c_str());
    if (!flist) {
        cerr << "ERROR: cannot open input list file -> " << input_list_filename << endl;
        return 1;
    }

    Int_t _numOfChannels = 0;
    Int_t _numOfFiles = 0;

    const Int_t MAX_N_channels = 4;
    const Int_t MAX_N_files = 1000;

    vector<vector<string>> list_filename(MAX_N_channels, vector<string>(MAX_N_files, ""));

    Int_t ch_idx = -1;
    Int_t file_idx = -1;
    Char_t Buffer[256];

    while (flist.getline(Buffer, sizeof(Buffer))) {
        if (Buffer[0] == '#') continue;
        istringstream strin(Buffer);
        string filename;
        if (!(strin >> ch_idx >> file_idx >> filename)) {
            cerr << "ERROR: invalid format in list file" << endl;
            return 1;
        }
        if (ch_idx >= MAX_N_channels || file_idx >= MAX_N_files) {
            cerr << "ERROR: overflow array for channels or files" << endl;
            return 1;
        }
        list_filename[ch_idx][file_idx] = filename;
        if (file_idx == 0) _numOfChannels++;
        if (ch_idx == 0) _numOfFiles++;
    }

    cout << "Channels: " << _numOfChannels << ", Files per channel: " << _numOfFiles << endl;

    // 出力ROOTファイルのオープン
    TFile* outFile = TFile::Open(output_root_filename.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        cerr << "ERROR: cannot open output ROOT file -> " << output_root_filename << endl;
        return 1;
    }

    // TTree の定義（生データのみ）
    TTree* tree = new TTree("tree", "Raw Waveform Tree");

    Int_t event;
    ULong64_t time_stamp;
    Int_t channel;
    UShort_t wave_raw[_DT5751Length];

    tree->Branch("event", &event, "event/I");
    tree->Branch("time_stamp", &time_stamp, "time_stamp/l");
    tree->Branch("channel", &channel, "channel/I");
    tree->Branch("wave_raw", wave_raw, Form("wave_raw[%d]/s", _DT5751Length));

    // イベントフィルターリストのロード、またはTQファイルからの自動カット抽出
    set<Int_t> target_events;
    bool use_event_filter = false;

    if (!target_mode.empty()) {
        if (target_mode == "fastn" || target_mode == "slown") {
            // coinファイルのパスを自動マッピング
            string basename = input_list_filename;
            size_t wave_pos = basename.find("wave");
            if (wave_pos != string::npos) {
                basename.replace(wave_pos, 4, "tq");
            }
            size_t last_slash_bin = basename.find_last_of("/\\");
            if (last_slash_bin != string::npos) {
                basename = basename.substr(last_slash_bin + 1);
            }
            size_t dot_pos = basename.find_last_of(".");
            if (dot_pos != string::npos) {
                basename = basename.substr(0, dot_pos);
            }

            // 候補ファイル名の作成 (tq や coincidence の挿入パターンの多様性に対応)
            vector<string> candidate_basenames;
            candidate_basenames.push_back(basename);
            if (basename.find("tq") == string::npos) {
                size_t v_pos = basename.find("_v");
                if (v_pos != string::npos) {
                    string inserted = basename;
                    inserted.insert(v_pos, "_tq");
                    candidate_basenames.push_back(inserted);
                }
                candidate_basenames.push_back(basename + "_tq");
            }

            vector<string> coin_paths_to_try;
            for (const auto& cb : candidate_basenames) {
                string filename_coin = cb + "_coincidence.root";
                coin_paths_to_try.push_back("root/" + filename_coin);
                if (!list_dir.empty()) {
                    coin_paths_to_try.push_back(list_dir + filename_coin);
                    coin_paths_to_try.push_back(list_dir + "../root/" + filename_coin);
                    coin_paths_to_try.push_back(list_dir + "root/" + filename_coin);
                }
            }

            TFile* f_coin = nullptr;
            string resolved_coin_path = "";
            for (const auto& p : coin_paths_to_try) {
                f_coin = TFile::Open(p.c_str(), "READ");
                if (f_coin && !f_coin->IsZombie()) {
                    resolved_coin_path = p;
                    break;
                }
                if (f_coin) {
                    delete f_coin;
                    f_coin = nullptr;
                }
            }

            if (f_coin) {
                TTree* t_coin = (TTree*)f_coin->Get("tree");
                if (t_coin) {
                    Int_t ev;
                    Double_t delta_T_us = 0.0;
                    if (target_mode == "fastn") {
                        t_coin->SetBranchAddress("fast_event", &ev);
                    } else {
                        t_coin->SetBranchAddress("slow_event", &ev);
                        t_coin->SetBranchAddress("delta_T_us", &delta_T_us);
                    }
                    Long64_t ents = t_coin->GetEntries();
                    for (Long64_t k = 0; k < ents; ++k) {
                        t_coin->GetEntry(k);
                        if (target_mode == "slown" && delta_T_us > 500.0) {
                            continue;
                        }
                        target_events.insert(ev);
                    }
                    use_event_filter = true;
                    cout << "Loaded " << target_events.size() << " target events for mode '" << target_mode 
                         << "' from auto-mapped coincidence file: " << resolved_coin_path << endl;
                } else {
                    cerr << "WARNING: Cannot find TTree 'tree' in mapped file -> " << resolved_coin_path << endl;
                }
                f_coin->Close();
                delete f_coin;
            } else {
                cerr << "WARNING: Cannot open auto-mapped coincidence ROOT file in any search locations. Decoding all events." << endl;
            }
        } else {
            // TQファイルのパスを自動マッピング
            string basename = input_list_filename;
            size_t wave_pos = basename.find("wave");
            if (wave_pos != string::npos) {
                basename.replace(wave_pos, 4, "tq");
            }
            size_t last_slash_bin = basename.find_last_of("/\\");
            if (last_slash_bin != string::npos) {
                basename = basename.substr(last_slash_bin + 1);
            }
            size_t dot_pos = basename.find_last_of(".");
            if (dot_pos != string::npos) {
                basename = basename.substr(0, dot_pos);
            }

            // 候補ファイル名の作成
            vector<string> candidate_basenames;
            candidate_basenames.push_back(basename);
            if (basename.find("tq") == string::npos) {
                size_t v_pos = basename.find("_v");
                if (v_pos != string::npos) {
                    string inserted = basename;
                    inserted.insert(v_pos, "_tq");
                    candidate_basenames.push_back(inserted);
                }
                candidate_basenames.push_back(basename + "_tq");
            }

            vector<string> tq_paths_to_try;
            for (const auto& cb : candidate_basenames) {
                string filename_tq = cb + ".root";
                tq_paths_to_try.push_back("root/" + filename_tq);
                if (!list_dir.empty()) {
                    tq_paths_to_try.push_back(list_dir + filename_tq);
                    tq_paths_to_try.push_back(list_dir + "../root/" + filename_tq);
                    tq_paths_to_try.push_back(list_dir + "root/" + filename_tq);
                }
            }

            TFile* f_tq = nullptr;
            string resolved_tq_path = "";
            for (const auto& p : tq_paths_to_try) {
                f_tq = TFile::Open(p.c_str(), "READ");
                if (f_tq && !f_tq->IsZombie()) {
                    resolved_tq_path = p;
                    break;
                }
                if (f_tq) {
                    delete f_tq;
                    f_tq = nullptr;
                }
            }

            if (f_tq) {
                TTree* t_tq = (TTree*)f_tq->Get("tree");
                if (t_tq) {
                    Int_t ev;
                    Double_t t0, t1;
                    t_tq->SetBranchAddress("event", &ev);
                    t_tq->SetBranchAddress("T0", &t0);
                    t_tq->SetBranchAddress("T1", &t1);
                    
                    Long64_t ents = t_tq->GetEntries();
                    for (Long64_t k = 0; k < ents; ++k) {
                        t_tq->GetEntry(k);
                        Double_t diff = t1 - t0;
                        if (target_mode == "gamma") {
                            if (diff >= 40.0 && diff < 60.0) {
                                target_events.insert(ev);
                            }
                        } else if (target_mode == "n") {
                            if (diff >= 60.0 && diff <= 120.0) {
                                target_events.insert(ev);
                            }
                        }
                    }
                    use_event_filter = true;
                    cout << "Loaded " << target_events.size() << " target events for mode '" << target_mode 
                         << "' (T1-T0 range) from auto-mapped TQ file: " << resolved_tq_path << endl;
                } else {
                    cerr << "WARNING: Cannot find TTree 'tree' in mapped TQ file -> " << resolved_tq_path << endl;
                }
                f_tq->Close();
                delete f_tq;
            } else {
                cerr << "WARNING: Cannot open auto-mapped TQ ROOT file in any search locations. Decoding all events." << endl;
            }
        }
    } else if (!event_list_txt.empty()) {
        ifstream fev(event_list_txt.c_str());
        if (fev) {
            Int_t ev;
            while (fev >> ev) {
                target_events.insert(ev);
            }
            use_event_filter = true;
            cout << "Loaded " << target_events.size() << " target events from filter list: " << event_list_txt << endl;
        } else {
            cerr << "WARNING: cannot open event filter list -> " << event_list_txt << ". Decoding all events." << endl;
        }
    }

    Char_t* buf = new Char_t[_DT5751DataSize];

    for (Int_t j = 0; j < _numOfFiles; j++) {
        vector<ifstream*> ifs(_numOfChannels, nullptr);

        // Open input binary files for this subset
        for (Int_t i = 0; i < _numOfChannels; i++) {
            ifs[i] = new ifstream();
            string bin_path = list_filename[i][j];
            
            // ファイル名部分（basename）のみを抽出
            string basename = bin_path;
            size_t last_slash_bin = bin_path.find_last_of("/\\");
            if (last_slash_bin != string::npos) {
                basename = bin_path.substr(last_slash_bin + 1);
            }
            
            // 探索するパス候補のリストを作成
            vector<string> paths_to_try;
            paths_to_try.push_back(bin_path); // 1. リスト記載のそのままのパス
            if (!list_dir.empty()) {
                paths_to_try.push_back(list_dir + bin_path);             // 2. リストと同じフォルダ直下の記載通りのパス
                paths_to_try.push_back(list_dir + "rawdata/" + bin_path); // 3. リストと同じフォルダの rawdata/ 配下の記載通りのパス
                paths_to_try.push_back(list_dir + "rawdata/" + basename); // 4. 【賢い吸収】リストと同じフォルダの rawdata/ 配下の純粋なファイル名
                paths_to_try.push_back(list_dir + basename);             // 5. 【賢い吸収】リストと同じフォルダの直下の純粋なファイル名
            }
            
            Bool_t opened = false;
            for (const auto& p : paths_to_try) {
                ifs[i]->open(p.c_str(), ios::binary);
                if (ifs[i]->is_open()) {
                    bin_path = p;
                    opened = true;
                    break;
                }
                ifs[i]->clear(); // エラーステートをクリアして次の試行へ
            }

            if (!opened) {
                // エラーで終了する代わりに、警告メッセージを出して is_open() = false のまま続行
                cout << "WARNING: File not found in any searched locations. Skipping channel " << i << ": " << basename << endl;
            }
        }

        cout << "Processing file subset " << j + 1 << "/" << _numOfFiles << "..." << endl;

        while (true) {
            Bool_t eof_flag = false;
            Int_t active_files = 0;
            
            for (Int_t i = 0; i < _numOfChannels; i++) {
                // オープンできていないファイル（チャンネル）はスキップ
                if (!ifs[i] || !ifs[i]->is_open()) {
                    continue;
                }
                active_files++;

                ifs[i]->read(buf, _DT5751DataSize);
                if (ifs[i]->eof()) {
                    eof_flag = true;
                    break;
                }

                DT5751WFdata* wf = (DT5751WFdata*)(buf);
                event = wf->header[4];
                time_stamp = (ULong64_t)(0x7fffffff & wf->header[5]);
                channel = i;

                // 生波形のコピー
                for (Int_t k = 0; k < _DT5751Length; k++) {
                    wave_raw[k] = wf->waveform[k];
                }

                // ターゲットリストに入っていないイベントはスキップ
                if (use_event_filter && target_events.find(event) == target_events.end()) {
                    continue;
                }

                // TTreeへの格納
                tree->Fill();
            }

            if (active_files == 0 || eof_flag) {
                if (eof_flag) {
                    // 各ファイルが終端に達したか確認
                    for (Int_t i = 0; i < _numOfChannels; i++) {
                        if (ifs[i] && ifs[i]->is_open() && !ifs[i]->eof()) {
                            ifs[i]->peek(); 
                            if (!ifs[i]->eof()) {
                                cerr << "WARNING: file size mismatch in file subset index " << j << endl;
                            }
                        }
                    }
                }
                break; // 次のファイルサブセットへ
            }
        }

        // Close files of this subset
        for (Int_t i = 0; i < _numOfChannels; i++) {
            if (ifs[i]) {
                if (ifs[i]->is_open()) {
                    ifs[i]->close();
                }
                delete ifs[i];
            }
        }
    }

    // ROOTファイルへの書き込みと終了
    outFile->Write();
    outFile->Close();
    delete outFile;
    delete[] buf;

    cout << "Finished exporting raw waveforms. ROOT file saved to: " << output_root_filename << endl;
    return 0;
}
