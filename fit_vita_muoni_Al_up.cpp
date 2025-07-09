#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMatrixDSym.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include <TPaveText.h>
#include <iomanip> 
#include <fstream>
#include <vector>
#include <string>

using namespace std;

// Doppio esponenziale + fondo costante
Double_t DoppioEsponenzialeFondoRetta(Double_t *x, Double_t *par) {
    return par[0] * exp(-x[0] / par[1]) + par[2] * exp(-x[0] / par[3]) + par[4];
}

void fit_vita_muoni_Al_up() {
    // --- Parametri base ---
    string file_name = "misura_vita_media_Al_25_04_14_dati_mu_alluminio_up.txt";
    const int nBins = 50;
    const double xMin = 1.4, xMax = 10.0;

    // --- Carica dati ---
    ifstream file(file_name);
    vector<double> data;
    double value;
    while (file >> value) data.push_back(value);
    file.close();

    // --- Canvas e istogramma ---
    TCanvas* c = new TCanvas("c_up", "Fit Muone Up in Alluminio", 1200, 1000);
    TH1F* h = new TH1F("Up Decays in Al", "", nBins, xMin, xMax);
    for (double x : data) h->Fill(x);
    h->Scale(1.0 / h->Integral());
    h->SetLineColor(kBlue);
    h->SetFillColorAlpha(kBlue, 0.3);
    h->SetLineWidth(1);
    h->GetXaxis()->SetTitle("#Deltat [#mu s]");
    h->GetYaxis()->SetTitle("Normalized Counts");
    h->Draw("HIST");

    // --- Funzione di fit ---
    TF1* fit = new TF1("fit_up", DoppioEsponenzialeFondoRetta, xMin, xMax, 5);
    fit->SetParameters(0.4, 2.2, 0.2, 0.864, 0.01);  // A, tau1, B, tau2, C
    fit->SetParNames("A", "#tau_{PVT}", "B", "#tau_{Al}", "C");
    fit->SetLineColor(kBlack);
    fit->SetLineWidth(2);

    // --- Vincoli ---
    fit->SetParLimits(0, 0, 100);      // A
    fit->SetParLimits(1, 2.1, 2.3);    // tau_PVT
    fit->SetParLimits(2, 0, 100);      // B
    fit->SetParLimits(3, 0.75, 0.9);   // tau_Al
    fit->SetParLimits(4, 0, 0.1);      // C

    // --- Fit con salvataggio del risultato ---
    gStyle->SetOptFit(1);
    TFitResultPtr fitResult = h->Fit(fit, "RS");  // "R" = range, "S" = save
    fit->Draw("SAME");

    // --- Etichette testo ---
    TLatex latex;
    latex.SetTextSize(0.035);
    latex.DrawLatexNDC(0.15, 0.83, Form("#tau_{PVT} = %.3f #pm %.3f #mus", fit->GetParameter(1), fit->GetParError(1)));
    latex.DrawLatexNDC(0.15, 0.77, Form("#tau_{Al} = %.3f #pm %.3f #mus", fit->GetParameter(3), fit->GetParError(3)));

    // --- Salvataggio grafico ---
    c->SaveAs("fit_muoni_Al_up.pdf");

    // === MATRICE DI CORRELAZIONE ===
    const int nPar = fit->GetNpar();
    TMatrixDSym corrMatrix = fitResult->GetCorrelationMatrix();

    // --- Salva in CSV ---
    ofstream csv("correlation_matrix_Al_up.csv");
    csv << std::fixed << std::setprecision(6);

    // Intestazione (nomi dei parametri)
    for (int j = 0; j < nPar; ++j) {
        csv << fit->GetParName(j);
        if (j < nPar - 1) csv << ",";
        else csv << "\n";
    }

    // Riga per riga della matrice
    for (int i = 0; i < nPar; ++i) {
        for (int j = 0; j < nPar; ++j) {
            csv << corrMatrix(i, j);
            if (j < nPar - 1) csv << ",";
        }
        csv << "\n";
    }

    csv.close();
}
