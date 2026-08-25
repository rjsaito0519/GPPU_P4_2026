#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>

using namespace std;

static const int _DT5751DataSize = 2090;
static const int _DT5751Length = 1029;

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
    string output_root_filename = "";

    // 引数の賢い判定
    if (argc == 3) {
        string arg2 = argv[2];
        // 第2引数の拡張子が .txt の場合はイベントリストとして扱う
        if (arg2.size() > 4 && arg2.substr(arg2.size() - 4) == ".txt") {
            event_list_txt = arg2;
        } else {
            output_root_filename = arg2;
        }
    } else if (argc > 3) {
        event_list_txt = argv[2];
        output_root_filename = argv[3];
    }

    // サフィックスの自動抽出 (例: data/Cf252_tq_merge_n_events.txt -> suffix = "_n")
    string suffix = "";
    if (!event_list_txt.empty()) {
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

    // 出力先ROOTファイル名の自動生成 (サフィックスを付与)
    if (output_root_filename.empty()) {
        output_root_filename = input_list_filename;
        size_t last_dot = output_root_filename.find_last_of(".");
        if (last_dot != string::npos) {
            output_root_filename = output_root_filename.substr(0, last_dot);
        }
        output_root_filename += suffix + ".root"; // 例: Cf252_wave_01_n.root などのサフィックス付き
    }

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

    int _numOfChannels = 0;
    int _numOfFiles = 0;

    const int MAX_N_channels = 4;
    const int MAX_N_files = 1000;

    vector<vector<string>> list_filename(MAX_N_channels, vector<string>(MAX_N_files, ""));

    int ch_idx = -1;
    int file_idx = -1;
    char Buffer[256];

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
    TTree* tree = new TTree("wave_tree", "Raw Waveform Tree");

    Int_t event;
    ULong64_t time_stamp;
    Int_t channel;
    UShort_t wave_raw[_DT5751Length];

    tree->Branch("event", &event, "event/I");
    tree->Branch("time_stamp", &time_stamp, "time_stamp/l");
    tree->Branch("channel", &channel, "channel/I");
    tree->Branch("wave_raw", wave_raw, Form("wave_raw[%d]/s", _DT5751Length));

    // イベントフィルターリストのロード
    set<Int_t> target_events;
    bool use_event_filter = false;
    if (!event_list_txt.empty()) {
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

    char* buf = new char[_DT5751DataSize];

    for (int j = 0; j < _numOfFiles; j++) {
        vector<ifstream*> ifs(_numOfChannels, nullptr);

        // Open input binary files for this subset
        for (int i = 0; i < _numOfChannels; i++) {
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
            
            bool opened = false;
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
            bool eof_flag = false;
            int active_files = 0;
            
            for (int i = 0; i < _numOfChannels; i++) {
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
                for (int k = 0; k < _DT5751Length; k++) {
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
                    for (int i = 0; i < _numOfChannels; i++) {
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
        for (int i = 0; i < _numOfChannels; i++) {
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
