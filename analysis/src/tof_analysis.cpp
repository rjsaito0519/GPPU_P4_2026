#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TF1.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TLine.h>

// フィット結果およびフィット関数のポインタを保持する構造体
struct FitParams {
    Double_t mean_g;
    Double_t mean_g_err;
    Double_t sigma_g;
    Double_t sigma_g_err;
    Double_t mean_n;
    Double_t mean_n_err;
    Double_t sigma_n;
    Double_t sigma_n_err;
    Double_t p2p;
    Double_t p2p_err;
    TF1* func_g;
    TF1* func_n;
};

// ピークフィッティングを行うヘルパー関数
FitParams fit_tof_peaks(TH1D* h, const std::string& suffix) {
    // --- 1. ガンマ線ピークのフィット (60 ns 未満) ---
    // Prefit 1回目 (広域: 40.0 ~ 58.0)
    TF1* f_pre1_g = new TF1(("f_pre1_g_" + suffix).c_str(), "gaus(0) + pol1(3)", 40.0, 58.0);
    f_pre1_g->SetParameters(h->GetMaximum() * 0.8, 50.0, 3.0, 10.0, 0.0);
    h->Fit(f_pre1_g, "R Q N");
    Double_t m1_g = f_pre1_g->GetParameter(1);
    Double_t s1_g = f_pre1_g->GetParameter(2);
    delete f_pre1_g;

    // Prefit 2回目 (m1 +- 3.0 * s1)
    Double_t p2_min_g = m1_g - 3.0 * s1_g;
    Double_t p2_max_g = m1_g + 3.0 * s1_g;
    TF1* f_pre2_g = new TF1(("f_pre2_g_" + suffix).c_str(), "gaus(0) + pol1(3)", p2_min_g, p2_max_g);
    f_pre2_g->SetParameters(h->GetMaximum() * 0.8, m1_g, s1_g, 10.0, 0.0);
    h->Fit(f_pre2_g, "R Q N");
    Double_t m2_g = f_pre2_g->GetParameter(1);
    Double_t s2_g = f_pre2_g->GetParameter(2);
    delete f_pre2_g;

    // 本フィット (m2 +- 1.5 * s2)
    Double_t fit_min_g = m2_g - 1.5 * s2_g;
    Double_t fit_max_g = m2_g + 1.5 * s2_g;
    TF1* f_g = new TF1(("f_g_" + suffix).c_str(), "gaus(0) + pol1(3)", fit_min_g, fit_max_g);
    f_g->SetParameters(h->GetMaximum() * 0.8, m2_g, s2_g, 10.0, 0.0);
    f_g->SetLineColor(kRed);
    f_g->SetLineWidth(3);
    h->Fit(f_g, "R Q N");

    // --- 2. 中性子ピークのフィット (60 ns 超) ---
    // 中性子は裾野が広いため prefit は1回のみとし、ガンマ線ピーク(60ns以下)を巻き込まないように左限を60.0nsに制限
    TF1* f_pre_n = new TF1(("f_pre_n_" + suffix).c_str(), "gaus(0) + pol1(3)", 62.0, 110.0);
    f_pre_n->SetParameters(h->GetMaximum() * 0.8, 80.0, 10.0, 20.0, 0.0);
    h->Fit(f_pre_n, "R Q N");
    Double_t m_n = f_pre_n->GetParameter(1);
    Double_t s_n = f_pre_n->GetParameter(2);
    delete f_pre_n;

    // 本フィット (m_n +- 1.5 * s_n)
    Double_t fit_min_n = std::max(60.0, m_n - 1.5 * s_n); // 左限を 60.0 ns にクリップ
    Double_t fit_max_n = m_n + 1.5 * s_n;
    TF1* f_n = new TF1(("f_n_" + suffix).c_str(), "gaus(0) + pol1(3)", fit_min_n, fit_max_n);
    f_n->SetParameters(h->GetMaximum() * 0.8, m_n, s_n, 20.0, 0.0);
    f_n->SetLineColor(kRed);
    f_n->SetLineWidth(3);
    h->Fit(f_n, "R Q N");

    FitParams p;
    p.mean_g = f_g->GetParameter(1);
    p.mean_g_err = f_g->GetParError(1);
    p.sigma_g = f_g->GetParameter(2);
    p.sigma_g_err = f_g->GetParError(2);
    p.mean_n = f_n->GetParameter(1);
    p.mean_n_err = f_n->GetParError(1);
    p.sigma_n = f_n->GetParameter(2);
    p.sigma_n_err = f_n->GetParError(2);
    p.p2p = p.mean_n - p.mean_g;
    p.p2p_err = std::sqrt(p.mean_g_err * p.mean_g_err + p.mean_n_err * p.mean_n_err);
    p.func_g = f_g;
    p.func_n = f_n;

    return p;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_root_path> [output_pdf_dir] [flight_path_length_m]" << std::endl;
        return 1;
    }

    std::string input_path = argv[1];
    std::string output_dir = "pdf"; // デフォルトの出力ディレクトリ
    Double_t flight_path = 1.0;     // デフォルトの飛行距離 (1.0 m)

    if (argc > 2) {
        output_dir = argv[2];
    }
    if (argc > 3) {
        try {
            flight_path = std::stod(argv[3]);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Invalid flight_path_length_m. Using default value 1.0 m." << std::endl;
            flight_path = 1.0;
        }
    }

    // 出力ディレクトリを自動作成 (なければ作成)
    gSystem->mkdir(output_dir.c_str(), true);

    TFile* file = TFile::Open(input_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open " << input_path << std::endl;
        return 1;
    }

    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: Cannot find TTree 'tree' in input file" << std::endl;
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
    std::string pdf_path = output_dir + "/" + base_name + "_tof_analysis.pdf";

    // 1〜4ページ用キャンバス (元のアスペクト比 800 x 600)
    TCanvas* c = new TCanvas("c", "TOF Analysis", 800, 600);

    // PDF の書き込み開始 (オープン)
    c->Print((pdf_path + "[").c_str());

    // -------------------------------------------------------------
    // 1. スルーイング補正用パラメータ p0 を自動決定するフィッティング
    // -------------------------------------------------------------
    std::cout << "Fitting slewing curve automatically..." << std::endl;
    TH2D* h2_fit = new TH2D("h2_fit", "Slewing Fit (T1-T0 in 30-60);Q0;T0", 100, 0, 50, 240, 280, 400);
    tree->Draw("T0:Q0>>h2_fit", "(T1-T0) < 60 && 30 < (T1-T0)", "goff");
    TProfile* prof = h2_fit->ProfileX("prof");
    TF1* f_slew = new TF1("f_slew", "[0]/sqrt(x) + [1]", 1.5, 50);
    f_slew->SetParameters(30.0, 350.0);

    prof->Fit(f_slew, "R Q");

    Double_t p0 = f_slew->GetParameter(0);
    Double_t p1 = f_slew->GetParameter(1);
    std::cout << "Determined Slewing Parameter p0: " << p0 << " (Offset p1: " << p1 << ")" << std::endl;

    h2_fit->SetTitle(Form("Slewing Fit (T1-T0 in 30-60) [p0 = %.2f];Q0;T0", p0));
    h2_fit->Draw("colz");
    prof->SetLineColor(kBlack);
    prof->SetLineWidth(2);
    prof->Draw("same");
    f_slew->SetLineColor(kRed);
    f_slew->SetLineWidth(3);
    f_slew->Draw("same");
    c->Print(pdf_path.c_str());

    // -------------------------------------------------------------
    // 2. T0 vs Q0 2Dヒストグラム (補正前)
    // -------------------------------------------------------------
    TH2D* h2_before = new TH2D("h2_before", "T0 vs Q0 (Before Correction) {Q1 > 0.0};Q0;T0", 150, 0, 50, 240, 280, 400);
    tree->Draw("T0:Q0>>h2_before", "Q1 > 0.0", "colz");
    c->Print(pdf_path.c_str());

    // -------------------------------------------------------------
    // 3. T0_corr vs Q0 2Dヒストグラム (その場で数式を用いて補正)
    // -------------------------------------------------------------
    TH2D* h2_after = new TH2D("h2_after", Form("T0_corr vs Q0 (After Correction) [p0=%.2f] {Q1 > 0.0};Q0;T0_corr", p0), 150, 0, 50, 240, 280, 400);
    tree->Draw(Form("T0 - (%f / sqrt(Q0)) : Q0 >> h2_after", p0), "Q1 > 0.0 && Q0 > 0.0", "colz");
    c->Print(pdf_path.c_str());

    // -------------------------------------------------------------
    // 4. TOF (T1 - T0) vs (T1 - T0_corr) の1D比較
    // -------------------------------------------------------------
    TH1D* h_tof_raw = new TH1D("h_tof_raw", "TOF Comparison {Q0 >= 4.0 && Q1 > 0.0};TOF [ns];Entries", 220, 20, 130);
    TH1D* h_tof_corr = new TH1D("h_tof_corr", "TOF Comparison {Q0 >= 4.0 && Q1 > 0.0};TOF [ns];Entries", 220, 20, 130);

    tree->Draw("T1 - T0>>h_tof_raw", "Q0 >= 4.0 && Q1 > 0.0", "goff");
    tree->Draw(Form("T1 - (T0 - (%f / sqrt(Q0)))>>h_tof_corr", p0), "Q0 >= 4.0 && Q1 > 0.0", "goff");

    h_tof_raw->SetLineColor(kBlue);
    h_tof_raw->SetLineWidth(2);
    h_tof_corr->SetLineColor(kRed);
    h_tof_corr->SetLineWidth(2);

    Double_t max_val = std::max(h_tof_raw->GetMaximum(), h_tof_corr->GetMaximum()) * 1.15;
    h_tof_raw->SetMaximum(max_val);

    h_tof_raw->Draw("hist");
    h_tof_corr->Draw("hist same");

    TLegend* leg = new TLegend(0.68, 0.76, 0.88, 0.88);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_tof_raw, "Before", "l");
    leg->AddEntry(h_tof_corr, "After", "l");
    leg->Draw();

    c->Print(pdf_path.c_str());

    // -------------------------------------------------------------
    // 5. Before 分布の個別フィッティング & パラメータ表示 (5ページ目)
    // -------------------------------------------------------------
    std::cout << "Fitting peaks for 'Before' TOF..." << std::endl;
    // 5〜6ページ用キャンバス (横長 1100 x 600)
    TCanvas* c_fit = new TCanvas("c_fit", "TOF Peak Fit", 1100, 600);
    c_fit->SetRightMargin(0.32); // 右マージンを広げて枠外にテキストを収める

    h_tof_raw->SetTitle("Before Correction: TOF Peak Fit;TOF [ns];Entries");
    h_tof_raw->SetLineColor(kBlack);
    h_tof_raw->Draw("hist");
    
    FitParams params_before = fit_tof_peaks(h_tof_raw, "before");
    params_before.func_g->Draw("same");
    params_before.func_n->Draw("same");

    // ピーク中心を示す縦点線 (赤色、スタイル2)
    TLine* line_g_before = new TLine(params_before.mean_g, 0, params_before.mean_g, h_tof_raw->GetMaximum() * 1.05);
    line_g_before->SetLineColor(kRed);
    line_g_before->SetLineStyle(2);
    line_g_before->SetLineWidth(2);
    line_g_before->Draw("same");

    TLine* line_n_before = new TLine(params_before.mean_n, 0, params_before.mean_n, h_tof_raw->GetMaximum() * 1.05);
    line_n_before->SetLineColor(kRed);
    line_n_before->SetLineStyle(2);
    line_n_before->SetLineWidth(2);
    line_n_before->Draw("same");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.032);
    latex.SetTextColor(kBlack); // すべて黒文字に統一
    latex.DrawLatex(0.70, 0.70, "#gamma peak:");
    latex.DrawLatex(0.70, 0.65, Form("  Mean : %.2f #pm %.2f", params_before.mean_g, params_before.mean_g_err));
    latex.DrawLatex(0.70, 0.60, Form("  Sigma: %.2f #pm %.2f", params_before.sigma_g, params_before.sigma_g_err));
    
    latex.DrawLatex(0.70, 0.50, "Neutron peak:");
    latex.DrawLatex(0.70, 0.45, Form("  Mean : %.2f #pm %.2f", params_before.mean_n, params_before.mean_n_err));
    latex.DrawLatex(0.70, 0.40, Form("  Sigma: %.2f #pm %.2f", params_before.sigma_n, params_before.sigma_n_err));
    
    latex.DrawLatex(0.70, 0.30, "Interval:");
    latex.DrawLatex(0.70, 0.25, Form("  P2P  : %.2f #pm %.2f", params_before.p2p, params_before.p2p_err));

    c_fit->Print(pdf_path.c_str());

    // -------------------------------------------------------------
    // 6. After 分布の個別フィッティング & パラメータ表示 (6ページ目)
    // -------------------------------------------------------------
    std::cout << "Fitting peaks for 'After' TOF..." << std::endl;
    c_fit->Clear();
    h_tof_corr->SetTitle("After Correction: TOF Peak Fit;TOF [ns];Entries");
    h_tof_corr->SetLineColor(kBlack);
    h_tof_corr->Draw("hist");
    
    FitParams params_after = fit_tof_peaks(h_tof_corr, "after");
    params_after.func_g->Draw("same");
    params_after.func_n->Draw("same");

    // ピーク中心を示す縦点線 (赤色、スタイル2)
    TLine* line_g_after = new TLine(params_after.mean_g, 0, params_after.mean_g, h_tof_corr->GetMaximum() * 1.05);
    line_g_after->SetLineColor(kRed);
    line_g_after->SetLineStyle(2);
    line_g_after->SetLineWidth(2);
    line_g_after->Draw("same");

    TLine* line_n_after = new TLine(params_after.mean_n, 0, params_after.mean_n, h_tof_corr->GetMaximum() * 1.05);
    line_n_after->SetLineColor(kRed);
    line_n_after->SetLineStyle(2);
    line_n_after->SetLineWidth(2);
    line_n_after->Draw("same");

    latex.SetTextColor(kBlack); // すべて黒文字に統一
    latex.DrawLatex(0.70, 0.70, "#gamma peak:");
    latex.DrawLatex(0.70, 0.65, Form("  Mean : %.2f #pm %.2f", params_after.mean_g, params_after.mean_g_err));
    latex.DrawLatex(0.70, 0.60, Form("  Sigma: %.2f #pm %.2f", params_after.sigma_g, params_after.sigma_g_err));
    
    latex.DrawLatex(0.70, 0.50, "Neutron peak:");
    latex.DrawLatex(0.70, 0.45, Form("  Mean : %.2f #pm %.2f", params_after.mean_n, params_after.mean_n_err));
    latex.DrawLatex(0.70, 0.40, Form("  Sigma: %.2f #pm %.2f", params_after.sigma_n, params_after.sigma_n_err));
    
    latex.DrawLatex(0.70, 0.30, "Interval:");
    latex.DrawLatex(0.70, 0.25, Form("  P2P  : %.2f #pm %.2f", params_after.p2p, params_after.p2p_err));

    c_fit->Print(pdf_path.c_str());

    // -------------------------------------------------------------
    // 7. 中性子運動エネルギー分布への変換 (7ページ目) [After 補正後のみ]
    // -------------------------------------------------------------
    std::cout << "Calculating and plotting Neutron Kinetic Energy (using After Correction)..." << std::endl;
    // 7ページ目も通常アスペクト比 800 x 600 に戻して見やすくします
    c->Clear();
    c->SetRightMargin(0.10); // 右マージンを元に戻す

    // 物理定数
    Double_t c_light = 0.299792458; // 光速 [m/ns]
    Double_t m_n = 939.565;         // 中性子の質量 [MeV]

    // ガンマ線ピーク位置と真のガンマ線TOF
    Double_t t0_gamma = params_after.mean_g;
    Double_t t_gamma = flight_path / c_light; // ガンマ線の真の飛行時間 [ns]

    // 横軸 0 ~ 15 MeV のヒストグラムを作成 (bin数: 75, 1bin=0.2MeV)
    TH1D* h_energy = new TH1D("h_energy", Form("Neutron Kinetic Energy Spectrum (L = %.2f m) {Q0 >= 4.0};Energy [MeV];Entries", flight_path), 75, 0, 15);
    h_energy->SetLineColor(kBlack);
    h_energy->SetLineWidth(2);

    // TTree::Draw 用の相対論的エネルギー変換式の構築
    // 1. T1 - T0_corr を求める数式:
    std::string t0_corr_str = Form("(T0 - (%f / sqrt(Q0)))", p0);
    std::string t_diff_str = "T1 - " + t0_corr_str;

    // 2. 真のTOF を求める数式: t_TOF = t_diff - t0_gamma + t_gamma
    std::string t_tof_str = Form("(%s) - %f + %f", t_diff_str.c_str(), t0_gamma, t_gamma);

    // 3. beta を求める数式: beta = L / (c_light * t_TOF)
    std::string beta_str = Form("%f / (%f * (%s))", flight_path, c_light, t_tof_str.c_str());

    // 4. Relative Kinetic Energy: E_k = (1 / sqrt(1 - beta^2) - 1) * m_n
    //    かつ、物理的イベント (0 < beta < 1) かつ (t_TOF > 0) のみを選択するフィルターをCut条件に含める
    std::string energy_formula = Form("(1.0 / sqrt(1.0 - (%s)*(%s)) - 1.0) * %f", beta_str.c_str(), beta_str.c_str(), m_n);
    
    // Cut条件: Q0 >= 4.0 && Q1 > 0.0 かつ 物理イベント (beta < 1.0 かつ t_TOF > 0.0)
    std::string cuts = Form("Q0 >= 4.0 && Q1 > 0.0 && (%s) < 1.0 && (%s) > 0.0", beta_str.c_str(), t_tof_str.c_str());

    tree->Draw((energy_formula + ">>h_energy").c_str(), cuts.c_str(), "hist");

    // Y軸の見切れを防ぐため、最大値を自動取得した最大値の 1.15倍に設定
    h_energy->SetMaximum(h_energy->GetMaximum() * 1.15);

    c->Print(pdf_path.c_str());

    // PDF の書き込み終了
    c->Print((pdf_path + "]").c_str());

    std::cout << "\nTOF Analysis completed successfully." << std::endl;
    std::cout << "Results saved to: " << pdf_path << std::endl;

    file->Close();
    delete params_before.func_g;
    delete params_before.func_n;
    delete params_after.func_g;
    delete params_after.func_n;
    delete file;
    delete c;
    delete c_fit;

    return 0;
}
