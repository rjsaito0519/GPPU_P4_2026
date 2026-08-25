#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <TFile.h>
#include <TTree.h>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_root_file> <output_txt_file> [tree_name]" << endl;
        return 1;
    }

    string input_root = argv[1];
    string output_txt = argv[2];
    string tree_name = "coincidence_tree";
    if (argc > 3) {
        tree_name = argv[3];
    }

    TFile* file = TFile::Open(input_root.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "ERROR: cannot open input ROOT file -> " << input_root << endl;
        return 1;
    }

    TTree* tree = (TTree*)file->Get(tree_name.c_str());
    if (!tree) {
        cerr << "ERROR: cannot find TTree '" << tree_name << "' in input file" << endl;
        file->Close();
        return 1;
    }

    Int_t event;
    tree->SetBranchAddress("event", &event);

    set<Int_t> unique_events;
    Long64_t n_entries = tree->GetEntries();
    
    cout << "Scanning TTree '" << tree_name << "' to collect event numbers..." << endl;
    for (Long64_t i = 0; i < n_entries; ++i) {
        tree->GetEntry(i);
        unique_events.insert(event);
    }

    ofstream ofs(output_txt.c_str());
    if (!ofs) {
        cerr << "ERROR: cannot open output text file -> " << output_txt << endl;
        file->Close();
        return 1;
    }

    for (const auto& ev : unique_events) {
        ofs << ev << "\n";
    }

    cout << "Successfully exported " << unique_events.size() << " unique events to: " << output_txt << endl;

    ofs.close();
    file->Close();
    delete file;
    return 0;
}
