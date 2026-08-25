#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

using namespace std;

static const int _DT5751DataSize = 2090;
static const int _DT5751Length = 1029;

typedef struct {
    unsigned int header[8];
    unsigned short waveform[_DT5751Length];
} DT5751WFdata;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_list_filename> <output_root_filename>" << endl;
        return 1;
    }

    string input_list_filename = argv[1];
    string output_root_filename = argv[2];

    ifstream flist(input_list_filename.c_str());
    if (!flist) {
        cerr << "ERROR: cannot open input list file -> " << input_list_filename << endl;
        return 1;
    }

    // パラメータ定義 (wave2tq.ccに準拠)
    const double threshold = 10.0; // ADC value
    const int n_baseline_length = 200;
    const int n_softtrigger_enable_length = 800;
    const int n_length_pre_softtrigger = 10;
    const int n_length_post_softtrigger = 100;
    const double impedance = 50.0; // ohm (nominal)
    const double dt = 1.0; // nsec (sampling interval)

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

    // TTree の定義
    TTree* tree = new TTree("wave_tree", "Waveform Tree");

    Int_t event;
    ULong64_t time_stamp;
    Int_t channel;
    UShort_t wave_raw[_DT5751Length];
    Float_t wave[_DT5751Length];
    Double_t baseline;
    Double_t charge;
    Double_t time;
    Int_t softtrigger;

    tree->Branch("event", &event, "event/I");
    tree->Branch("time_stamp", &time_stamp, "time_stamp/l");
    tree->Branch("channel", &channel, "channel/I");
    tree->Branch("wave_raw", wave_raw, Form("wave_raw[%d]/s", _DT5751Length));
    tree->Branch("wave", wave, Form("wave[%d]/F", _DT5751Length));
    tree->Branch("baseline", &baseline, "baseline/D");
    tree->Branch("charge", &charge, "charge/D");
    tree->Branch("time", &time, "time/D");
    tree->Branch("softtrigger", &softtrigger, "softtrigger/I");

    char* buf = new char[_DT5751DataSize];

    for (int j = 0; j < _numOfFiles; j++) {
        vector<ifstream*> ifs(_numOfChannels, nullptr);

        // Open input binary files for this subset
        for (int i = 0; i < _numOfChannels; i++) {
            ifs[i] = new ifstream(list_filename[i][j].c_str(), ios::binary);
            if (!ifs[i] || !ifs[i]->is_open()) {
                cerr << "ERROR: cannot open input binary file -> " << list_filename[i][j] << endl;
                return 1;
            }
        }

        cout << "Processing file subset " << j + 1 << "/" << _numOfFiles << "..." << endl;

        while (true) {
            bool eof_flag = false;
            
            for (int i = 0; i < _numOfChannels; i++) {
                ifs[i]->read(buf, _DT5751DataSize);
                if (ifs[i]->eof()) {
                    eof_flag = true;
                    break;
                }

                DT5751WFdata* wf = (DT5751WFdata*)(buf);
                event = wf->header[4];
                time_stamp = (ULong64_t)(0x7fffffff & wf->header[5]);
                channel = i;

                // 1. Raw waveform copy
                for (int k = 0; k < _DT5751Length; k++) {
                    wave_raw[k] = wf->waveform[k];
                }

                // 2. Baseline estimation (wave2tq.cc algorithm)
                double sum_wave = 0.0;
                int n_wave = 0;
                for (int k = 0; k < n_baseline_length; k++) {
                    sum_wave += (double)wave_raw[k];
                    n_wave++;
                }
                double baseline_rough = sum_wave / double(n_wave);

                sum_wave = 0.0;
                n_wave = 0;
                for (int k = 0; k < n_baseline_length; k++) {
                    if (abs((double)wave_raw[k] - baseline_rough) < threshold) {
                        sum_wave += (double)wave_raw[k];
                        n_wave++;
                    }
                }
                baseline = (n_wave > 0) ? (sum_wave / double(n_wave)) : baseline_rough;

                // 3. Baseline subtraction (invert polarity)
                for (int k = 0; k < _DT5751Length; k++) {
                    wave[k] = (Float_t)(baseline - (double)wave_raw[k]);
                }

                // 4. Threshold trigger search & Time determination
                charge = 0.0;
                softtrigger = 0;
                int k_threshold = 0;

                for (int k = 0; k < n_softtrigger_enable_length; k++) {
                    if (wave[k] >= threshold) {
                        softtrigger = 1;
                        k_threshold = k;
                        break;
                    }
                }

                // 5. Charge integration
                if (softtrigger == 1) {
                    int k_window_start = k_threshold - n_length_pre_softtrigger;
                    int k_window_end = k_threshold + n_length_post_softtrigger;

                    for (int k = max(0, k_window_start); k < min(_DT5751Length, k_window_end); k++) {
                        charge += ((double)wave[k] / impedance * dt); // pC
                    }
                }

                // 6. Threshold constant fraction / linear interpolation timing
                time = 0.0;
                if (k_threshold > 0) {
                    double time_1 = double(k_threshold - 1);
                    double time_2 = double(k_threshold);
                    double wave_1 = (double)wave[k_threshold - 1];
                    double wave_2 = (double)wave[k_threshold];

                    double slope = (wave_1 - wave_2) / (time_1 - time_2);
                    double offset = (time_1 * wave_2 - time_2 * wave_1) / (time_1 - time_2);

                    if (slope != 0.0) {
                        time = (threshold - offset) / slope;
                    }
                }

                // Fill event to tree
                tree->Fill();
            }

            if (eof_flag) {
                // Confirm all files in the subset reached EOF
                for (int i = 0; i < _numOfChannels; i++) {
                    if (!ifs[i]->eof()) {
                        ifs[i]->peek(); 
                        if (!ifs[i]->eof()) {
                            cerr << "WARNING: file size mismatch in file subset index " << j << endl;
                        }
                    }
                }
                break; // End of subset
            }
        }

        // Close files of this subset
        for (int i = 0; i < _numOfChannels; i++) {
            ifs[i]->close();
            delete ifs[i];
        }
    }

    // ROOTファイルへの書き込みと終了
    outFile->Write();
    outFile->Close();
    delete outFile;
    delete[] buf;

    cout << "Finished exporting waveforms. ROOT file saved to: " << output_root_filename << endl;
    return 0;
}
