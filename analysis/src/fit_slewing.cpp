#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <iostream>

void fit_slewing(const char* filename = "data/Cf252_tq_01.root") {
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open " << filename << std::endl;
        return;
    }
    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: Cannot find tree" << std::endl;
        file->Close();
        return;
    }

    // T0 vs Q0 の2次元ヒストグラムを作成 (ユーザー指定のT1-T0 ゲートを適用: 30 ~ 60 ns)
    // 見切れを防ぐため、Y軸の範囲を 280 ~ 400 に広げます (240 bins)
    TH2D* h2 = new TH2D("h2", "T0 vs Q0 {(T1-T0) < 60 && 30 < (T1-T0)};Q0;T0", 100, 0, 50, 240, 280, 400);
    tree->Draw("T0:Q0>>h2", "(T1-T0) < 60 && 30 < (T1-T0)", "goff");

    // X軸（Q0）スライスごとの平均値プロファイルを求める
    TProfile* prof = h2->ProfileX("prof");

    // スルーイング補正用フィット関数: f(q) = p0 / sqrt(q) + p1
    // Q0が非常に小さい領域（1.5以下）はノイズや波形整形閾値によるカットオフで不安定なため、フィット範囲を 1.5 ~ 50 に設定
    TF1* f_slew = new TF1("f_slew", "[0]/sqrt(x) + [1]", 1.5, 50);
    f_slew->SetParameters(30.0, 350.0); // 初期値の設定

    std::cout << "Fitting profile..." << std::endl;
    prof->Fit(f_slew, "R");

    Double_t p0 = f_slew->GetParameter(0);
    Double_t p1 = f_slew->GetParameter(1);

    std::cout << "\n==============================================" << std::endl;
    std::cout << " Slewing Correction Parameters Determined!" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Suggested Formula: T0_corr = T0 - (" << p0 << " / sqrt(Q0))" << std::endl;
    std::cout << "Parameter p0 (amplitude)  : " << p0 << std::endl;
    std::cout << "Parameter p1 (offset constant): " << p1 << std::endl;
    std::cout << "==============================================\n" << std::endl;

    TCanvas* c1 = new TCanvas("c1", "Slewing Fit", 800, 600);
    h2->Draw("colz");
    prof->SetLineColor(kBlack);
    prof->SetLineWidth(2);
    prof->Draw("same");
    f_slew->SetLineColor(kRed);
    f_slew->SetLineWidth(3);
    f_slew->Draw("same");

    c1->SaveAs("data/slewing_fit.png");
    std::cout << "Plot saved to data/slewing_fit.png" << std::endl;

    file->Close();
    delete file;
}

// 実行バイナリ用メイン関数
int main(int argc, char** argv) {
    const char* filename = "data/Cf252_tq_01.root";
    if (argc > 1) {
        filename = argv[1];
    }
    fit_slewing(filename);
    return 0;
}
