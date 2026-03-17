// templateFitRoo_BDT.C
#include "TFile.h"
#include "TH1.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace RooFit;


void TemplateFitterDb(const char* dataFile="data.root", const char* dataHist="hist",
                        const char* mc_b_File="template_b.root", const char* mc_b_Hist="hist_b",
                        const char* mc_c_File="template_c.root", const char* mc_c_Hist="hist_c",
                        const char* mc_lf_File="template_lf.root", const char* mc_lf_Hist="hist_lf",
                        TString ptbin = "5-10") {

  TFile *fData = TFile::Open(dataFile);
  TH1 *hData = (TH1*)fData->Get(dataHist);

  TFile *fb = TFile::Open(mc_b_File);
  TH1 *hB = (TH1*)fb->Get(mc_b_Hist);

  TFile *fc = TFile::Open(mc_c_File);
  TH1 *hC = (TH1*)fc->Get(mc_c_Hist);

  TFile *fl = TFile::Open(mc_lf_File);
  TH1 *hL = (TH1*)fl->Get(mc_lf_Hist);

  hData->Scale(1./hData->Integral());
  hB->Scale(1./hB->Integral());
  hC->Scale(1./hC->Integral());
  hL->Scale(1./hL->Integral());

  if(!hData || !hB || !hC || !hL) {
    printf("Missing histograms\n");
    return;
  }
    const double eps = 1e-20;

  auto protect = [&](TH1* h){
    for(int ix=1; ix<=h->GetNbinsX(); ++ix){
      if(h->GetBinContent(ix) <= 0)
        h->SetBinContent(ix, eps);
    }
  };

  protect(hB);
  protect(hC);
  protect(hL);

  double xmin = hData->GetXaxis()->GetXmin();
  double xmax = hData->GetXaxis()->GetXmax();
  RooRealVar x("x","D_bscore", xmin, xmax);

  RooDataHist dataHistR("dataHistR","data binned", RooArgList(x), Import(*hData));
  RooDataHist bHistR("bHistR","b template", RooArgList(x), Import(*hB));
  RooDataHist cHistR("cHistR","c template", RooArgList(x), Import(*hC));
  RooDataHist lHistR("lHistR","lf template", RooArgList(x), Import(*hL));

  RooArgSet  obsSet(x);

  RooHistPdf pdfB("pdfB","pdfB", obsSet, bHistR);
  RooHistPdf pdfC("pdfC","pdfC", obsSet, cHistR);
  RooHistPdf pdfL("pdfL","pdfL", obsSet, lHistR);

  // double Ndata = hData->Integral();
  // RooRealVar yB("yB","yield B", 0.01, 0.0, 1);
  // RooRealVar yC("yC","yield C", 0.03, 0.0, 1);
  // RooRealVar yL("yL","yield L", 0.96, 0.0, 1);

  RooRealVar fB("fB","fraction B", 0.01, 0.0, 1.0);
  RooRealVar fC("fC","fraction C", 0.1, 0.0, 1.0);

  // fL is implicit: 1 - fB - fC
  RooFormulaVar fL("fL","1 - fB - fC", RooArgList(fB,fC));

  // if (ptbin == "5-10") fC.setVal(0.02);
  // if (ptbin == "10-20") fC.setVal(0.1);
  // if (ptbin == "20-40") fC.setVal(0.2);
  // if (ptbin == "40-70") fC.setVal(0.35);
  // if (ptbin == "70-120") fC.setVal(0.5);
  // if (ptbin == "120-250") fC.setVal(0.6);

  // fC.setConstant(kTRUE);

  RooAddPdf model("model","model",
                  RooArgList(pdfB,pdfC,pdfL),
                  RooArgList(fB,fC));
  //RooAddPdf model("model","model", RooArgList(pdfB,pdfC,pdfL), RooArgList(yB,yC,yL));

  RooFitResult *res = model.fitTo(dataHistR, Extended(kTRUE), Save(kTRUE), Offset(kTRUE), SumW2Error(kTRUE));

  res->Print();

    RooPlot *frame = x.frame();
    RooCurve *c_data = (RooCurve*) dataHistR.plotOn(frame, Name("data"));
    model.plotOn(frame, LineColor(kRed), LineWidth(2), Name("model_fit"));
    model.plotOn(frame, Components("pdfB"), LineColor(kBlue), LineStyle(kDashed), Name("b_plot"));
    model.plotOn(frame, Components("pdfC"), LineColor(kGreen+2), LineStyle(kDashed), Name("c_plot"));
    model.plotOn(frame, Components("pdfL"), LineColor(kMagenta), LineStyle(kDashed), Name("lf_plot"));

    RooHist *hist_data = (RooHist *) frame->findObject("data");
    RooCurve *c_fit = (RooCurve*) frame->findObject("model_fit");
    RooCurve *c_b = (RooCurve*) frame->findObject("b_plot");
    RooCurve *c_c = (RooCurve*) frame->findObject("c_plot");
    RooCurve *c_lf = (RooCurve*) frame->findObject("lf_plot");

    TLine *line = new TLine(-10.,1., +10., 1.);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->SetLineStyle(7);

    TH1* hModel = (TH1*)hData->Clone("hModel");

    for (int ix = 1; ix <= hData->GetNbinsX(); ++ix){
      double xcenter = hData->GetXaxis()->GetBinCenter(ix);
      x.setVal(xcenter);
      double val = model.getVal(obsSet);
      hModel->SetBinContent(ix, val);
    }
    hModel->Scale(1./hModel->Integral());
    TH1 *hRatio = (TH1*)hModel->Clone("hRatio");
    
    hRatio->Divide(hData, hModel, 1., 1., "B");

  TCanvas *c = new TCanvas("c","template fit",800,900);
    gStyle->SetOptStat(0000);
  float topPadFrac = 0.70;
  TPad *pad1 = new TPad("pad1","Main Pad",0,1-topPadFrac,1,1);
  TPad *pad2 = new TPad("pad2","Ratio Pad",0,0,1,1-topPadFrac);

  pad1->Draw();
  pad2->Draw();

  float leftMargin   = 0.14;
  float rightMargin  = 0.04;
  float topMargin    = 0.05;
  float bottomMargin = 0.35;

  pad1->SetLeftMargin(leftMargin);
  pad1->SetRightMargin(rightMargin);
  pad1->SetTopMargin(topMargin);
  pad1->SetBottomMargin(0.001);
  pad1->SetLogy();
  pad1->SetTicks(1);

  pad2->SetLeftMargin(leftMargin);
  pad2->SetRightMargin(rightMargin);
  pad2->SetTopMargin(0.001);
  pad2->SetBottomMargin(bottomMargin);
  pad2->SetTicks(1);

  pad1->cd();
  frame->GetYaxis()->SetRangeUser(7e-6, 2e-1);
  frame->SetTitle("");
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->SetTitle("D_{b}");
  frame->GetYaxis()->SetTitle("Probability Density");
  frame->Draw();


TLegend *legend = new TLegend(0.7,0.45,0.98,0.65);
legend->SetBorderSize(0);
legend->SetFillStyle(0);
legend->SetTextFont(42);

legend->AddEntry(c_b, "b jets", "L");
legend->AddEntry(c_c, "c jets", "L");
legend->AddEntry(c_lf, "lf jets", "L");
legend->AddEntry(hist_data, "Data", "LEP");
legend->AddEntry(c_fit, "Fit", "L");

legend->Draw();

TPaveText *expt = new TPaveText(0.5,0.7,0.95,0.89,"NDC");
expt->AddText("pp, #sqrt{#it{s}} = 13.6 TeV");
expt->AddText("Anti-#it{k}_{T}, R = 0.4, |#eta| < 0.5");
expt->AddText(Form("%s GeV/#it{c}",ptbin.Data()));

expt->SetBorderSize(0);
expt->SetFillStyle(0);
expt->SetTextFont(42);

expt->Draw();

c->Update();

pad2->cd();
  hRatio->GetYaxis()->SetRangeUser(0.55, 1.45);
  hRatio->GetYaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetTitleSize(0.12);
  hRatio->GetYaxis()->SetTitleOffset(0.6);
  hRatio->GetXaxis()->SetLabelSize(0.12);
  hRatio->GetXaxis()->SetTitleSize(0.15);
  hRatio->GetXaxis()->SetTickSize(0.05);
  hRatio->GetYaxis()->SetNdivisions(7,5,0);
  hRatio->GetXaxis()->SetTitle("D_{b}");
  hRatio->GetYaxis()->SetTitle("Data/Fit");
  hRatio->SetMarkerStyle(kFullCircle);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetLineColor(kBlack);
  hRatio->Draw();
  line->Draw();

  // double total = yB.getVal() + yC.getVal() + yL.getVal();
  // printf("Fractions: B=%.4f +/- %.4f, C=%.4f +/- %.4f, L=%.4f +/- %.4f\n",yB.getVal()/total, yB.getError()/total,yC.getVal()/total, yC.getError()/total,yL.getVal()/total, yL.getError()/total);
  // printf("Abundance ratio = %.4f : %.4f : %.4f\n", yB.getVal()/yB.getVal(), yC.getVal()/yB.getVal(), yL.getVal()/yB.getVal());

    double fracb = fB.getVal();
  double fracc = fC.getVal();
  double fracl = 1.0 - fracb - fracc;

  printf("\nFractions:\n");
  printf("B = %.4f\n", fracb);
  printf("C = %.4f\n", fracc);
  printf("L = %.4f\n", fracl);

  Printf("Abundance Ratios : %f : %f : %f\n", fracb/fracb, fracc/fracb, fracl/fracb);


  c->SaveAs("Db_score_template_fit.pdf");

  TH1D *fractions = new TH1D("fractions", "", 3, 0, 3);
  fractions->SetBinContent(1, fL.getVal());
  fractions->SetBinError(1, TMath::Sqrt(TMath::Power(fB.getError(),2.) + TMath::Power(fC.getError(),2.)));
  
  fractions->SetBinContent(2, fC.getVal());
  fractions->SetBinError(2, fC.getError());

  fractions->SetBinContent(3, fB.getVal());
  fractions->SetBinError(3, fB.getError());

  // TH1D *fractions = new TH1D("fractions", "", 3, 0, 3);
  // fractions->SetBinContent(1, yL.getVal()/total);
  // fractions->SetBinError(1, yL.getError()/total);
  
  // fractions->SetBinContent(2, yC.getVal()/total);
  // fractions->SetBinError(2, yC.getError()/total);

  // fractions->SetBinContent(3, yB.getVal()/total);
  // fractions->SetBinError(3, yB.getError()/total);

  TFile *fitresult = new TFile("fitresultDb.root", "recreate");
  fitresult->cd();

  res->Write();
  c->Write();
  fitresult->Write();
  fractions->Write();
  hRatio->Write();
  fitresult->Close();
}
