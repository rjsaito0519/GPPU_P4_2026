#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TAxis.h>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_psd_root_file>" << endl;
        return 1;
    }

    string input_path = argv[1];

    TFile* file = TFile::Open(input_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "ERROR: cannot open input ROOT file -> " << input_path << endl;
        return 1;
    }

    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        cerr << "ERROR: cannot find TTree 'tree' in input file" << endl;
        file->Close();
        return 1;
    }

    TApplication app("app", &argc, argv);

    TCanvas* c = new TCanvas("c_fom", "PSD Fit and FOM Evaluation", 800, 600);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetGrid();

    // PSD のヒストグラムを作成 (0.0 から 0.6 の範囲)
    TH1D* h_psd = new TH1D("h_psd", "PSD Spectrum and Double-Gaussian Fit;PSD;Entries", 120, 0.0, 0.6);
    h_psd->SetLineColor(kBlack);
    h_psd->SetLineWidth(2);

    // ツリーからPSDをプロット
    tree->Draw("PSD >> h_psd", "PSD > 0.0");

    // 初期パラメータの自動推定
    // ガンマピークの初期推定 (0.08 ~ 0.16 の範囲の最大値)
    int bin_min_g = h_psd->FindBin(0.08);
    int bin_max_g = h_psd->FindBin(0.16);
    double max_g = 0;
    double mean_g_init = 0.12;
    for (int b = bin_min_g; b <= bin_max_g; ++b) {
        if (h_psd->GetBinContent(b) > max_g) {
            max_g = h_psd->GetBinContent(b);
            mean_g_init = h_psd->GetBinCenter(b);
        }
    }

    // 中性子ピークの初期推定 (0.20 ~ 0.35 の範囲の最大値)
    int bin_min_n = h_psd->FindBin(0.20);
    int bin_max_n = h_psd->FindBin(0.35);
    double max_n = 0;
    double mean_n_init = 0.25;
    for (int b = bin_min_n; b <= bin_max_n; ++b) {
        if (h_psd->GetBinContent(b) > max_n) {
            max_n = h_psd->GetBinContent(b);
            mean_n_init = h_psd->GetBinCenter(b);
        }
    }

    // ダブルガウスフィット関数の定義: gaus(0) + gaus(3)
    // パラメータ: [0]:Amp_g, [1]:Mean_g, [2]:Sigma_g, [3]:Amp_n, [4]:Mean_n, [5]:Sigma_n
    TF1* f_double = new TF1("f_double", "gaus(0) + gaus(3)", 0.02, 0.55);
    
    f_double->SetParameter(0, max_g);
    f_double->SetParameter(1, mean_g_init);
    f_double->SetParameter(2, 0.02); // ガンマ線の幅の初期値は狭め
    f_double->SetParLimits(2, 0.005, 0.06);

    f_double->SetParameter(3, max_n);
    f_double->SetParameter(4, mean_n_init);
    f_double->SetParameter(5, 0.05); // 中性子の幅の初期値は広め
    f_double->SetParLimits(5, 0.01, 0.12);

    f_double->SetLineColor(kRed);
    f_double->SetLineWidth(3);

    cout << "Fitting PSD spectrum..." << endl;
    h_psd->Fit(f_double, "R");

    // フィット結果の抽出
    double amp_g   = f_double->GetParameter(0);
    double mean_g  = f_double->GetParameter(1);
    double sigma_g = f_double->GetParameter(2);
    double amp_n   = f_double->GetParameter(3);
    double mean_n  = f_double->GetParameter(4);
    double sigma_n = f_double->GetParameter(5);

    // FWHM の計算
    double fwhm_g = 2.355 * sigma_g;
    double fwhm_n = 2.355 * sigma_n;

    // FOM (Figure of Merit) の計算
    double fom = 0.0;
    if ((fwhm_g + fwhm_n) > 0.0) {
        fom = abs(mean_n - mean_g) / (fwhm_g + fwhm_n);
    }

    cout << "\n=============================================" << endl;
    cout << "  PSD Double-Gaussian Fit Results" << endl;
    cout << "=============================================" << endl;
    cout << " Gamma Peak (left):" << endl;
    cout << "  - Mean (P_gamma):  " << mean_g << "  [initial: " << mean_g_init << "]" << endl;
    cout << "  - FWHM_gamma:      " << fwhm_g << " (Sigma: " << sigma_g << ")" << endl;
    cout << " Neutron Peak (right):" << endl;
    cout << "  - Mean (P_n):      " << mean_n << "  [initial: " << mean_n_init << "]" << endl;
    cout << "  - FWHM_neutron:    " << fwhm_n << " (Sigma: " << sigma_n << ")" << endl;
    cout << "---------------------------------------------" << endl;
    cout << "  >>> Figure of Merit (FOM): " << fom << " <<<" << endl;
    cout << "=============================================\n" << endl;

    // 個々のガウス成分を描画して視覚的に確認できるようにする
    TF1* f_g = new TF1("f_g", "gaus", 0.0, 0.6);
    f_g->SetParameters(amp_g, mean_g, sigma_g);
    f_g->SetLineColor(kBlue);
    f_g->SetLineStyle(2);
    f_g->Draw("same");

    TF1* f_n = new TF1("f_n", "gaus", 0.0, 0.6);
    f_n->SetParameters(amp_n, mean_n, sigma_n);
    f_n->SetLineColor(kGreen+2);
    f_n->SetLineStyle(2);
    f_n->Draw("same");

    c->Update();

    // 出力PDFフォルダの作成と保存
    string base_name = input_path;
    size_t last_slash = base_name.find_last_of("/\\");
    if (last_slash != string::npos) {
        base_name = base_name.substr(last_slash + 1);
    }
    size_t last_dot = base_name.find_last_of(".");
    if (last_dot != string::npos) {
        base_name = base_name.substr(0, last_dot);
    }
    string out_pdf = "pdf/" + base_name + "_fom_fit.pdf";
    gSystem->mkdir("pdf", true);
    c->Print(out_pdf.c_str());
    cout << "Saved fit plot to: " << out_pdf << endl;

    app.Run();

    file->Close();
    delete file;
    return 0;
}
