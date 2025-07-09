/**
 * @author      Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        ex1.C
 * @created     Mon May 08, 2023 10:32:12 CEST
 */

#include <iostream>
#include "TSystem.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TTimer.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMatrixT.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/DataOptions.h"
#include "Fit/Chi2FCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "Math/Functor.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"


/*
const int    data_nbin = 100; 
const double data_xmin = 0; 
const double data_xmax = 100; 
const double data_exposure = 1000; 
const double data_bin_width = (data_xmax - data_xmin) / data_nbin; 
*/
//typedef std::vector<ParameterDef_t> parset;
/*
double fc_bkg(const double* x, const double* par) {
  const double norm = 1.0 / (par[1]-par[1]*exp(-data_xmax/par[1]) );
  return par[0] * norm * exp(-x[0] / par[1]) * data_bin_width; 
}

double fc_sig(const double* x, const double* par) {
  return par[0]*TMath::Gaus(x[0], par[1], par[2], true) * data_bin_width; 
}
*/

double fc_exp(const double* x, const double* par) {
  return par[0] * exp(-x[0] / par[1]); 
}

/*
TH1D* build_toy_dataset(TF1* fModel, const double* par);
*/
struct PullTerm_t {
  int fIdx; 
  double fValue; 
  double fErr; 

  PullTerm_t(int i, double v, double err) : fIdx(i), fValue(v), fErr(err) {}; 
}; 

struct GlobalLkl_t{
  GlobalLkl_t(ROOT::Math::IMultiGenFunction &f1, std::vector<PullTerm_t>& pulls) 
    : fFCN(&f1), fPenalties(pulls) {}

  double operator()(const double *par) const
  {
    double penalty = eval_penalty(par); 
    return (*fFCN)(par) + penalty;
  }

  double eval_penalty(const double* params) const {
    double penalty = 0; 
    for (const auto &pull : fPenalties) {
      double par = params[pull.fIdx]; 
      penalty += 0.5*TMath::Sq( (par - pull.fValue) / pull.fErr); 
    }

    return penalty;
  }

  private: 
  const ROOT::Math::IMultiGenFunction* fFCN; 
  std::vector<PullTerm_t> fPenalties;
};

void ex2b(const bool debug = true) 
{

  string file_name = "misura_vita_media_Al_25_04_14_dati_mu_alluminio_down.txt";

  const int nBins = 50;
  const double xMin = 1.4, xMax = 10.0;

  // --- Carica dati ---
  ifstream file(file_name);
  vector<double> data;
  double value;
  while (file >> value) data.push_back(value);
  file.close();

  // --- Canvas e istogramma ---
 // TCanvas* c = new TCanvas("c_up", "Fit Muone Down in Alluminio", 1200, 1000);
  TH1F* h = new TH1F("Down Decays in Al", "", nBins, xMin, xMax);
  for (double x : data) h->Fill(x);
 // h->Scale(1.0 / h->Integral());
  h->SetLineColor(kRed);
  h->SetFillColorAlpha(kRed, 0.3);
  h->SetLineWidth(1);
  h->GetXaxis()->SetTitle("#Deltat [#mu s]");
  h->GetYaxis()->SetTitle("Normalized Counts");
  h->Draw("HIST");

  const double par_init[5] = {1000, 2.2, 600, 0.864, 50}; 
  const TString par_name[5] = {"N_PVT", "tau_PVT", "N_Al", "tau_Al", "Bkg"};

  // define model for background and compute normalization in the analysis domain
  TF1* modelDecayPVT = new TF1("modelDecayPVT", fc_exp, xMin, xMax, 2); 
  modelDecayPVT->SetParameters(par_init[0], par_init[1]); 

  TF1* modelDecayAl = new TF1("modelDecayAl", fc_exp, xMin, xMax, 2); 
  modelDecayAl->SetParameters(par_init[2], par_init[3]); 

  // define model for signal and compute normalization in the analysis domain
  TF1* modelBkg = new TF1("modelBgk", "pol0", xMin, xMax); 
  modelBkg->SetParameters(par_init[4]); 


  // define data model (sig + bkg)
  TF1* modelData = new TF1("modelData", 
      [&](double* x, double *p){
      double y = fc_exp(x, p) + fc_exp(x, &p[2]) + p[4];  
      return y;},
      xMin, xMax, 5); 
  modelData->SetParName(0, "N_{PVT}"); 
  modelData->SetParName(1, "#tau_{PVT}"); 
  modelData->SetParName(2, "N_{Al}"); 
  modelData->SetParName(3, "#tau_{Al}"); 
  modelData->SetParName(4, "Bkg"); 

  modelData->SetParameters( par_init ); 
  const int n_par = modelData->GetNpar(); 
/*
  parset params(n_par); 
  parset truth(n_par);
  for (int i=0; i<n_par; i++) {
    ParameterDef_t pdef; 
    pdef.fName  = par_name[i]; 
    pdef.fTitle = modelData->GetParName(i); 
    pdef.fVal = 0.; 
    params.at(i) = pdef; 
    pdef.fVal = par_init[i]; 
    truth.at(i) = pdef;
  }
*/
  //const double par_init[5] = {1.53, 28, 0.5, 30, 7};


  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange fit_range(1.4, 10); 
  ROOT::Math::WrappedMultiTF1 model(*modelData, n_par); 

  ROOT::Fit::Fitter fitter; 
  fitter.Config().SetParamsSettings(5, par_init);
  fitter.Config().MinimizerOptions().SetErrorDef(0.5); 
  fitter.Config().MinimizerOptions().SetPrintLevel(1);
  fitter.Config().SetMinimizer("Minuit2", "Migrad");

  for (int n=0; n<n_par; n++) {
    fitter.Config().ParSettings(n).SetName(modelData->GetParName(n)); 
  }
 
  // create containers to store fit results
  /*
  TFile* output_file = new TFile("ex2b_output.root", "recreate"); 
  TTree* tree = new TTree("best_fit", "Best fit results"); 
  for (auto &p: params) {
    tree->Branch(p.fName, &p.fVal); 
    tree->Branch(p.fName+"_err", &p.fStdDev); 
  }
*/
  TCanvas* cDebug = nullptr; 
  //TTimer* timer = nullptr;
  //if (debug) {
  //  timer = new TTimer("gSystem->ProcessEvents();", 500, false);
    cDebug = new TCanvas("cDebug", "debug", 0, 0, 1000, 600); 
  //  cDebug->Divide(2, 1); 
 // }

  //auto fit_toy = [&](TH1D* h, parset& params) {
    ROOT::Fit::BinData data_hist(opt, fit_range);
    ROOT::Fit::FillData(data_hist, h); 
    ROOT::Fit::PoissonLLFunction fcn(data_hist, model);

    std::vector<PullTerm_t> pulls; 
    pulls.push_back( PullTerm_t(1, 2.197, 0.05)); 

    GlobalLkl_t lkl(fcn, pulls); 

    //printf("hist counts = %g\n", data_hist.SumOfContent()); 

    int status = fitter.FitFCN(n_par, lkl, par_init, data_hist.Size(), false); 

    auto fit_result = fitter.Result();
/*
    for (int i=0; i<n_par; i++) {
      ParameterDef_t* pdef = &params[i]; 
      pdef->fVal = fit_result.Parameter(i); 
      pdef->fStdDev = fit_result.ParError(i); 
    }
*/
    /*if (debug) {
      timer->TurnOn(); 
      timer->Reset(); 
      cDebug->Clear("D"); 
      cDebug->cd(1); 
      */
      
      for (int i=0; i<n_par; i++) {
        modelData->SetParameter(i, fit_result.Parameter(i)   ); 
        modelData->SetParError (i, fit_result.ParError(i)    ); 
      }

      for(int i = 0; i < n_par; i++) {
        double Err_low, Err_up;
        fitter.GetMinimizer()->GetMinosError(i, Err_low, Err_up);
        printf("%s BestFit: %g (ErrFit: %g) - Minos: %g - %g\n", modelData->GetParName(i), modelData->GetParameter(i), modelData->GetParError(i), Err_low, Err_up);
      }

      gStyle->SetOptFit(1);
      h->GetListOfFunctions()->Add(modelData);
      h->Draw("hist"); 
    //  modelData->Draw("same"); 
/*
      cDebug->cd(2); 
      unsigned int npx = 100; 
      double xval[npx]; double yval[npx]; 
      fit_result.Scan(4, npx, xval, yval, 
          params.at(4).fVal - 3*params.at(4).fStdDev,
          params.at(4).fVal + 3*params.at(4).fStdDev);
      double ymin = *(std::min_element(yval, yval+npx)); 
      TGraph gProfile(npx); 
      TGraph gPull(npx); 
      for (int i=0; i<npx; i++) {
        gProfile.SetPoint(i, xval[i], yval[i] - ymin); 
        gPull.SetPoint(i, xval[i], 
            0.5*pow( (xval[i] - pulls.front().fValue) / pulls.front().fErr, 2) ); 
      }
      gProfile.SetLineWidth(2); gProfile.SetLineColor(kRed+1); 
      gPull.SetLineWidth(2); gPull.SetLineColor(kBlue+1); gPull.SetLineStyle(7);  
      gProfile.Draw("awl"); 
      gPull.Draw("l");
*/
/*
      cDebug->cd(); 
      cDebug->Modified(); 
      cDebug->Update(); 
      getchar( ); 
      timer->TurnOff();
    }
*/
  //  return fit_result;
  //};

/*
  for (int j=0; j<5e3; j++) {
    auto h = build_toy_dataset(modelData, par_true);

    auto fit_result = fit_toy(h, params); 

    tree->Fill(); 


    delete h; 
  }


  tree->Write(); 

  gStyle->SetPalette(kBlackBody); 
  TColor::InvertPalette(); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  TCorrPlot corr_plot(900, 900); 
  corr_plot.SetMargin(80); 
  for (const auto &p : params) {
    corr_plot.RegisterPar(p); 
  }

  corr_plot.LoadTree(tree); 

  corr_plot.Fill(); 

  auto cplot = corr_plot.MakePlot(); 

  corr_plot.PrintCorr(); 
  corr_plot.PrintGaussFit(); 

  cplot->Draw(); 

  //output_file->Close(); 
*/

  return;
}

/*
TH1D* build_toy_dataset(TF1* fModel, const double* par) {
  TH1D* h = new TH1D("h", "h", data_nbin, data_xmin, data_xmax); 

  fModel->SetParameters( par ); 

  const double n_events_exp = data_exposure*(
      fModel->GetParameter(0) + fModel->GetParameter(2)); 
  int n_sample = gRandom->Poisson( n_events_exp ); 

  for (int i=0; i<n_sample; i++) { h->Fill( fModel->GetRandom() ); }

  return h;
}
*/

