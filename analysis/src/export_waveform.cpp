#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
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
        cerr << "Usage: " << argv[0] << " <input_list_filename> [output_root_filename]" << endl;
        return 1;
    }

    string input_list_filename = argv[1];
    string output_root_filename = "";
    if (argc > 2) {
        output_root_filename = argv[2];
    } else {
        output_root_filename = input_list_filename;
        size_t last_dot = output_root_filename.find_last_of(".");
        if (last_dot != string::npos) {
            output_root_filename = output_root_filename.substr(0, last_dot);
        }
        output_root_filename += ".root";
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
