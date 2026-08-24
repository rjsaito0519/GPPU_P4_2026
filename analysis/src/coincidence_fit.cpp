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
    
    // PDF の書き込み開始 (オープン)
    c->Print((pdf_path + "[").c_str());

    // Q1 ゲート設定: 20までは5刻み、それ以降50までは10刻み
    std::vector<std::pair<Double_t, Double_t>> gates = {
        {0.0, 5.0}, {5.0, 10.0}, {10.0, 15.0}, {15.0, 20.0},
        {20.0, 30.0}, {30.0, 40.0}, {40.0, 50.0}
    };
    
    std::vector<Double_t> vec_q1;
    std::vector<Double_t> vec_q1_err;
    std::vector<Double_t> vec_tau;
    std::vector<Double_t> vec_tau_err;

    c->Clear();
    c->Divide(2, 2);
    int pad_idx = 1;

    for (const auto& gate : gates) {
        Double_t q_min = gate.first;
        Double_t q_max = gate.second;
        std::string gate_cut = Form("slow_Q1 >= %f && slow_Q1 < %f", q_min, q_max);
        
        // 統計チェック
        Long64_t entries = tree->GetEntries(gate_cut.c_str());
        if (entries < 30) {
            std::cout << "Q1 in [" << q_min << ", " << q_max << "]: skipped due to low statistics (" << entries << " entries)" << std::endl;
            continue;
        }

        c->cd(pad_idx);
        
        // 0 ~ 800 us の範囲で 100 ビン (1bin = 8.0 us)
        std::string hist_name = Form("h_q1_%d_%d", (int)q_min, (int)q_max);
        TH1D* h = new TH1D(hist_name.c_str(), Form("slow_Q1: %d to %d;#Delta t [#mus];Entries", (int)q_min, (int)q_max), 100, 0, 800);
        h->SetLineColor(kBlack);
        h->SetLineWidth(2);

        tree->Draw(Form("delta_T_us>>%s", hist_name.c_str()), gate_cut.c_str(), "goff");

        // -----------------------------------------------------------------
        // Y軸の最大値を調整: 0 ~ 30 us 付近の巨大ノイズを除外した30us以降の最大値の 1.25 倍に設定
        // (1bin = 8.0us なので、30us は 3.75bin目付近。安全のため5bin目(32us~)以降の最大値を探す)
        // -----------------------------------------------------------------
        Double_t local_max = 0.0;
        for (int bin = 5; bin <= h->GetNbinsX(); ++bin) {
            Double_t content = h->GetBinContent(bin);
            if (content > local_max) {
                local_max = content;
            }
        }
        h->SetMaximum(local_max * 1.25);

        // フィット関数: [0]*exp(-x/[1]) + [2] (指数関数 + 定数項)
        // フィット範囲: 30 ~ 780 us
        TF1* f_exp = new TF1(Form("f_exp_%s", hist_name.c_str()), "[0]*exp(-x/[1]) + [2]", 30.0, 780.0);
        
        // 500 ~ 800 us (70bin から 100bin) の平均値をバックグラウンド初期値とする
        Double_t bg_est = 0.0;
        int bg_bins = 0;
        for (int bin = 70; bin <= h->GetNbinsX(); ++bin) {
            bg_est += h->GetBinContent(bin);
            bg_bins++;
        }
        if (bg_bins > 0) bg_est /= bg_bins;

        // 時定数 tau の初期値は、実際のスケールに合わせた 120.0 us を設定
        f_exp->SetParameters(local_max - bg_est, 120.0, bg_est);
        f_exp->SetParLimits(0, 0.0, local_max * 2.0); // 振幅に安全制限を設定
        f_exp->SetParLimits(1, 1.0, 700.0);           // 時定数の上限を 700 us に設定
        f_exp->FixParameter(2, bg_est);               // 背景を固定してフィッティングを劇的に安定化

        // 2段階フィットに L オプション (ポアソン対数尤度フィット) を追加して誤差を安定化
        h->Fit(f_exp, "R Q N L");
        h->Fit(f_exp, "R Q L");

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
        if (tau > 2.0 && tau < 700.0 && tau_err < tau * 0.5) {
            vec_q1.push_back(q_min + (q_max - q_min) / 2.0);
            vec_q1_err.push_back((q_max - q_min) / 2.0);
            vec_tau.push_back(tau);
            vec_tau_err.push_back(tau_err);
        }

        pad_idx++;
        if (pad_idx > 4) {
            c->Print(pdf_path.c_str());
            c->Clear();
            c->Divide(2, 2);
            pad_idx = 1;
        }
    }

    if (pad_idx > 1) {
        // 未使用のパッドをクリア
        for (int p = pad_idx; p <= 4; ++p) {
            c->cd(p)->Clear();
        }
        c->Print(pdf_path.c_str());
    }

    // -------------------------------------------------------------
    // 合算フィット (Total Coincidence Fit) [X+1ページ目]
    // -------------------------------------------------------------
    std::cout << "Plotting Total (all Q1) Coincidence Fit..." << std::endl;
    c->Clear();
    c->SetRightMargin(0.10); // 余白をデフォルトに戻す

    TH1D* h_total = new TH1D("h_total", "Total Coincidence Decay Fit (slow_Q1: 0 to 50);#Delta t [#mus];Entries", 100, 0, 800);
    h_total->SetLineColor(kBlack);
    h_total->SetLineWidth(2);

    tree->Draw("delta_T_us>>h_total", "slow_Q1 >= 0 && slow_Q1 < 50", "goff");

    // Y軸の最大値を調整 (30us/5bin以降の最大値の 1.25倍)
    Double_t total_max = 0.0;
    for (int bin = 5; bin <= h_total->GetNbinsX(); ++bin) {
        Double_t content = h_total->GetBinContent(bin);
        if (content > total_max) {
            total_max = content;
        }
    }
    h_total->SetMaximum(total_max * 1.25);

    TF1* f_total = new TF1("f_total", "[0]*exp(-x/[1]) + [2]", 30.0, 780.0);
    Double_t bg_total = 0.0;
    int bg_total_bins = 0;
    for (int bin = 70; bin <= h_total->GetNbinsX(); ++bin) {
        bg_total += h_total->GetBinContent(bin);
        bg_total_bins++;
    }
    if (bg_total_bins > 0) bg_total /= bg_total_bins;

    f_total->SetParameters(total_max - bg_total, 120.0, bg_total);
    f_total->SetParLimits(0, 0.0, total_max * 2.0);
    f_total->SetParLimits(1, 1.0, 700.0);
    f_total->FixParameter(2, bg_total); // 背景を固定

    // 2段階フィット
    h_total->Fit(f_total, "R Q N L");
    h_total->Fit(f_total, "R Q L");

    h_total->Draw("hist");
    f_total->SetLineColor(kRed);
    f_total->SetLineWidth(3);
    f_total->Draw("same");

    TLatex latex_total;
    latex_total.SetNDC();
    latex_total.SetTextSize(0.045);
    latex_total.SetTextColor(kRed+2);
    latex_total.DrawLatex(0.45, 0.75, Form("#tau = %.1f #pm %.1f #mus", f_total->GetParameter(1), f_total->GetParError(1)));
    latex_total.DrawLatex(0.45, 0.68, Form("BG = %.1f #pm %.1f", f_total->GetParameter(2), f_total->GetParError(2)));

    c->Print(pdf_path.c_str());

    // -------------------------------------------------------------
    // 最終ページ: 時定数 tau の Q1 依存性グラフを描画
    // -------------------------------------------------------------
    std::cout << "Plotting summary graph..." << std::endl;
    c->Clear();

    if (!vec_q1.empty()) {
        TGraphErrors* gr = new TGraphErrors(vec_q1.size(), &vec_q1[0], &vec_tau[0], &vec_q1_err[0], &vec_tau_err[0]);
        gr->SetTitle("Decay Constant #tau vs slow_Q1;slow_Q1;Decay Constant #tau [#mus]");
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.5);
        gr->SetMarkerColor(kBlue);
        gr->SetLineColor(kBlue);
        gr->SetLineWidth(2);
        
        // 縦軸レンジの自動最適化 (得られた最大のtauに応じてマージンを調整、最低でも 50.0 us は確保)
        Double_t max_tau = *std::max_element(vec_tau.begin(), vec_tau.end());
        gr->GetYaxis()->SetRangeUser(0.0, std::max(50.0, max_tau * 1.25));
        gr->GetXaxis()->SetRangeUser(0.0, 50.0); // Q1の上限は50まで
        
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
