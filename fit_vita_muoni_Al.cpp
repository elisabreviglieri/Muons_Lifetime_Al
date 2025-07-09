#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include <fstream>
#include <vector>
#include <string>

using namespace std;

// Funzione singolo esponenziale + fondo costante
Double_t EsponenzialeFondoRetta(Double_t *x, Double_t *par) {
    return par[0] * exp(-x[0] / par[1]) + par[2];
}

// Funzione doppio esponenziale + fondo costante
Double_t DoppioEsponenzialeFondoRetta(Double_t *x, Double_t *par) {
    return par[0] * exp(-x[0] / par[1]) + par[2] * exp(-x[0] / par[3]) + par[4];
}

// Estrazione delle componenti esponenziali
vector<TF1*> extract_exp(const TF1* fit, double xmin, double xmax) {
    vector<TF1*> output;
    if (fit->GetNpar() == 3) {
        TF1* fexp = new TF1("exp", "[0]*exp(-x/[1])", xmin, xmax);
        fexp->SetParameters(fit->GetParameter(0), fit->GetParameter(1));
        fexp->SetLineStyle(2);
        fexp->SetLineColor(kGray+2);
        output.push_back(fexp);
    } else if (fit->GetNpar() == 5) {
        TF1* f1 = new TF1("exp1", "[0]*exp(-x/[1])", xmin, xmax);
        TF1* f2 = new TF1("exp2", "[0]*exp(-x/[1])", xmin, xmax);
        f1->SetParameters(fit->GetParameter(0), fit->GetParameter(1));
        f2->SetParameters(fit->GetParameter(2), fit->GetParameter(3));
        f1->SetLineStyle(2); f2->SetLineStyle(2);
        f1->SetLineColor(kBlack); f2->SetLineColor(kGray+1);
        output.push_back(f1);
        output.push_back(f2);
    }
    return output;
}

void fit_vita_muoni_Al() {
    vector<string> file_names = {
        "misura_vita_media_Al_25_04_14_dati_mu_alluminio_up.txt",
        "misura_vita_media_Al_25_04_14_dati_mu_alluminio_down.txt"
    };
    vector<string> labels = {"Up", "Down"};
    vector<Color_t> colors = {kRed, kBlue};

    const int nBins = 50;
    const double xMin = 1.4, xMax = 16.0;

    for (int i = 0; i < 2; ++i) {
        // Leggi i dati
        ifstream file(file_names[i]);
        vector<double> data;
        double value;
        while (file >> value) data.push_back(value);
        file.close();

        // Canvas
        TCanvas* c = new TCanvas(Form("c%d", i), Form("Fit %s", labels[i].c_str()), 1200, 1000);

        // Istogramma
        TH1F* h = new TH1F(Form("h%d", i), Form("Muon Lifetime %s", labels[i].c_str()), nBins, xMin, xMax);
        for (double x : data) h->Fill(x);
        h->Scale(1.0 / h->Integral());
        h->SetLineColor(colors[i]);
        h->SetMarkerColor(colors[i]);
        h->SetMarkerStyle(20);
        h->Draw();

        // Fit doppio esponenziale con parametri e limiti ottimizzati
        TF1* fit = new TF1(Form("fit%d", i), DoppioEsponenzialeFondoRetta, xMin, xMax, 5);
        fit->SetParameters(0.4, 2.2, 0.2, 0.864, 0.01); // A1, tau1, A2, tau2, fondo
        fit->SetParLimits(0, 0, 100);        // A1
        fit->SetParLimits(1, 1.5, 2.5);      // tau1 (muoni liberi)
        fit->SetParLimits(2, 0, 100);        // A2
        fit->SetParLimits(3, 0.7, 1.0);      // tau2 (muoni legati in Al)
        fit->SetParLimits(4, 0, 0.1);        // fondo costante

        h->Fit(fit, "R");
        gStyle->SetOptFit(1);

        // Legenda
        TLegend* leg = new TLegend(0.65, 0.75, 0.90, 0.90);
        leg->AddEntry(h, Form("Muone %s", labels[i].c_str()), "lep");
        leg->AddEntry(fit, "Fit complessivo", "l");
        leg->Draw();

        // Componenti esponenziali
        auto components = extract_exp(fit, xMin, xMax);
        for (auto& fexp : components) fexp->Draw("same");

        // Latex con tau1 e tau2
        double tau1 = fit->GetParameter(1);
        double err1 = fit->GetParError(1);
        double tau2 = fit->GetParameter(3);
        double err2 = fit->GetParError(3);
        TLatex latex;
        latex.SetTextSize(0.035);
        latex.DrawLatexNDC(0.15, 0.83, Form("#tau_{1} = %.3f #pm %.3f #mus", tau1, err1));
        latex.DrawLatexNDC(0.15, 0.77, Form("#tau_{2} = %.3f #pm %.3f #mus", tau2, err2));
    }
}
