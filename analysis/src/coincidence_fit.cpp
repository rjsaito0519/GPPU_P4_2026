#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLatex.h>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <coincidence_root_path> [output_pdf_dir]" << std::endl;
        return 1;
    }

    std::string input_path = argv[1];
    std::string output_dir = "pdf";
    if (argc > 2) {
        output_dir = argv[2];
    }

    // 出力ディレクトリを自動作成 (なければ作成)
    gSystem->mkdir(output_dir.c_str(), true);

    TFile* file = TFile::Open(input_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open " << input_path << std::endl;
        return 1;
    }

    TTree* tree = (TTree*)file->Get("coincidence_tree");
    if (!tree) {
        std::cerr << "Error: Cannot find TTree 'coincidence_tree' in input file" << std::endl;
        file->Close();
        return 1;
    }

    // スタイルの共通設定
    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");

    // 出力PDFパスの組み立て
    size_t last_slash = input_path.find_last_of("/\\");
    std::string base_name = (last_slash == std::string::npos) ? input_path : input_path.substr(last_slash + 1);
    size_t last_dot = base_name.find_last_of(".");
    if (last_dot != std::string::npos) {
        base_name = base_name.substr(0, last_dot);
    }
    std::string pdf_path = output_dir + "/" + base_name + "_fit_results.pdf";

    // 2x2 分割キャンバス (1200 x 900)
    TCanvas* c = new TCanvas("c", "Coincidence Fit", 1200, 900);
    c->Print((pdf_path + "[").c_str());

    // -------------------------------------------------------------
    // 1ページ目: 先発事象 (fast) の Time Window (fast_T1 - fast_T0) を表示
    // -------------------------------------------------------------
    std::cout << "Plotting Trigger Time Window..." << std::endl;
    c->Clear();
    c->SetRightMargin(0.10); // 余白を通常サイズに設定
    
    TH1D* h_window = new TH1D("h_window", "Trigger Time Window (fast_T1 - fast_T0);fast_T1 - fast_T0 [ns];Entries", 100, 50, 130);
    h_window->SetLineColor(kBlue);
    h_window->SetLineWidth(2);
    h_window->SetFillColor(kBlue);
    h_window->SetFillStyle(3002);
    
    tree->Draw("fast_T1 - fast_T0>>h_window", "", "hist");
    h_window->Draw("hist");
    
    c->Print(pdf_path.c_str());
    delete h_window;

    // -------------------------------------------------------------
    // 2ページ目以降: Q1 5刻みごとの時間差分布フィット
    // -------------------------------------------------------------
    const Double_t q1_step = 5.0;
    const Double_t q1_max_limit = 100.0;
    
    std::vector<Double_t> vec_q1;
    std::vector<Double_t> vec_q1_err;
    std::vector<Double_t> vec_tau;
    std::vector<Double_t> vec_tau_err;

    c->Clear();
    c->Divide(2, 2);
    int pad_idx = 1;

    for (Double_t q_min = 0.0; q_min < q1_max_limit; q_min += q1_step) {
        Double_t q_max = q_min + q1_step;
        std::string gate_cut = Form("slow_Q1 >= %f && slow_Q1 < %f", q_min, q_max);
        
        // 統計チェック
        Long64_t entries = tree->GetEntries(gate_cut.c_str());
        if (entries < 30) { // 統計下限を30に緩和
            std::cout << "Q1 in [" << q_min << ", " << q_max << "]: skipped due to low statistics (" << entries << " entries)" << std::endl;
            continue;
        }

        c->cd(pad_idx);
        
        // 0 ~ 500 us の範囲で 100 ビン (1bin = 5 us)
        std::string hist_name = Form("h_q1_%d_%d", (int)q_min, (int)q_max);
        TH1D* h = new TH1D(hist_name.c_str(), Form("slow_Q1: %d to %d;#Delta t [#mus];Entries", (int)q_min, (int)q_max), 100, 0, 500);
        h->SetLineColor(kBlack);
        h->SetLineWidth(2);

        tree->Draw(Form("delta_T_us>>%s", hist_name.c_str()), gate_cut.c_str(), "goff");

        // フィット関数: [0]*exp(-x/[1]) + [2] (指数関数 + 定数項)
        // フィット範囲: 20 ~ 450 us (即発ノイズを避けるため 20 us から開始し、500us近辺の端を避けるため 450 us まで)
        TF1* f_exp = new TF1(Form("f_exp_%s", hist_name.c_str()), "[0]*exp(-x/[1]) + [2]", 20.0, 450.0);
        
        Double_t max_val = h->GetMaximum();
        Double_t bg_est = h->GetBinContent(95); // 475 us 付近の値をバックグラウンド初期値とする
        f_exp->SetParameters(max_val - bg_est, 30.0, bg_est);
        f_exp->SetParLimits(1, 1.0, 300.0);

        h->Fit(f_exp, "R Q");

        Double_t tau = f_exp->GetParameter(1);
        Double_t tau_err = f_exp->GetParError(1);

        h->Draw("hist");
        f_exp->SetLineColor(kRed);
        f_exp->SetLineWidth(3);
        f_exp->Draw("same");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.045);
        latex.SetTextColor(kRed+2);
        latex.DrawLatex(0.45, 0.75, Form("#tau = %.1f #pm %.1f #mus", tau, tau_err));
        latex.DrawLatex(0.45, 0.68, Form("BG = %.1f #pm %.1f", f_exp->GetParameter(2), f_exp->GetParError(2)));

        // グラフ用データに追加 (フィット結果が妥当な場合のみ)
        if (tau > 2.0 && tau < 300.0 && tau_err < tau * 0.5) {
            vec_q1.push_back(q_min + q1_step / 2.0);
            vec_q1_err.push_back(q1_step / 2.0);
            vec_tau.push_back(tau);
            vec_tau_err.push_back(tau_err);
        }

        pad_idx++;
        if (pad_idx > 4) {
            c->Print(pdf_path.c_str());
            c->Clear("D");
            c->Divide(2, 2);
            pad_idx = 1;
        }
    }

    if (pad_idx > 1) {
        for (int p = pad_idx; p <= 4; ++p) {
            c->cd(p)->Clear();
        }
        c->Print(pdf_path.c_str());
    }

    // -------------------------------------------------------------
    // 最終ページ: 時定数 tau の Q1 依存性グラフを描画
    // -------------------------------------------------------------
    std::cout << "Plotting summary graph..." << std::endl;
    c->Clear();
    c->SetRightMargin(0.10);

    if (!vec_q1.empty()) {
        TGraphErrors* gr = new TGraphErrors(vec_q1.size(), &vec_q1[0], &vec_tau[0], &vec_q1_err[0], &vec_tau_err[0]);
        gr->SetTitle("Decay Constant #tau vs slow_Q1;slow_Q1;Decay Constant #tau [#mus]");
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.5);
        gr->SetMarkerColor(kBlue);
        gr->SetLineColor(kBlue);
        gr->SetLineWidth(2);
        
        gr->GetYaxis()->SetRangeUser(0.0, 150.0);
        gr->GetXaxis()->SetRangeUser(0.0, q1_max_limit);
        
        gr->Draw("AP");
        c->Print(pdf_path.c_str());
        delete gr;
    } else {
        std::cerr << "Warning: No valid fit points to plot on summary graph." << std::endl;
    }

    // PDF 終了
    c->Print((pdf_path + "]").c_str());

    std::cout << "Coincidence fit completed. Results saved to: " << pdf_path << std::endl;

    file->Close();
    delete file;
    delete c;

    return 0;
}
