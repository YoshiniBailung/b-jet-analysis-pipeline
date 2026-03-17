#include "TFile.h"
#include "TH2.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooArgSet.h"

using namespace RooFit;

void TemplateFitter2D(const char* dataFile="data_2d.root",
                      const char* dataHist="hist",
                      const char* mc_b_File="template_2d_b.root",
                      const char* mc_b_Hist="hist_b",
                      const char* mc_c_File="template_2d_c.root",
                      const char* mc_c_Hist="hist_c",
                      const char* mc_lf_File="template_2d_lf.root",
                      const char* mc_lf_Hist="hist_lf",
                      TString ptbin = "5-10") {


  TFile *fData = TFile::Open(dataFile);
  TFile *fb = TFile::Open(mc_b_File);
  TFile *fc = TFile::Open(mc_c_File);
  TFile *fl = TFile::Open(mc_lf_File);

  TH2 *hData = (TH2*)fData->Get(dataHist);
  TH2 *hB    = (TH2*)fb->Get(mc_b_Hist);
  TH2 *hC    = (TH2*)fc->Get(mc_c_Hist);
  TH2 *hL    = (TH2*)fl->Get(mc_lf_Hist);

  if(!hData || !hB || !hC || !hL) {
    printf("Missing histograms\n");
    return;
  }

  hData->Scale(1./hData->Integral());
  hB->Scale(1./hB->Integral());
  hC->Scale(1./hC->Integral());
  hL->Scale(1./hL->Integral());

  const double eps = 1e-20;

  auto protect = [&](TH2* h){
    for(int ix=1; ix<=h->GetNbinsX(); ++ix){
      for(int iy=1; iy<=h->GetNbinsY(); ++iy){
        if(h->GetBinContent(ix,iy) <= 0)
          h->SetBinContent(ix,iy, eps);
      }
    }
  };

  protect(hB);
  protect(hC);
  protect(hL);

  // --------------------------
  // Define observables
  // --------------------------

  double xmin = hData->GetXaxis()->GetXmin();
  double xmax = hData->GetXaxis()->GetXmax();
  double ymin = hData->GetYaxis()->GetXmin();
  double ymax = hData->GetYaxis()->GetXmax();

  RooRealVar x("x","BDT_bclfscore", xmin, xmax);
  RooRealVar y("y","BDT_bcscore", ymin, ymax);

  RooArgList obs(x,y);
  RooArgSet  obsSet(x,y);

  // --------------------------
  // Import histograms
  // --------------------------

  RooDataHist dataR("dataR","data", obs, Import(*hData));
  RooDataHist bR("bR","b template", obs, Import(*hB));
  RooDataHist cR("cR","c template", obs, Import(*hC));
  RooDataHist lR("lR","lf template", obs, Import(*hL));



  RooHistPdf pdfB("pdfB","pdfB", obsSet, bR);
  RooHistPdf pdfC("pdfC","pdfC", obsSet, cR);
  RooHistPdf pdfL("pdfL","pdfL", obsSet, lR);


  RooRealVar fB("fB","fraction B", 0.03, 0.0, 1.0);
  RooRealVar fC("fC","fraction C", 0.1, 0.0, 1.0);
  //RooRealVar fL("fL","fraction L", 0.67, 0.0, 1.0);
  //RooGaussian cConstraint("cConstraint","c constraint",fC, 0.3, 0.1);

  // fL is implicit: 1 - fB - fC
  RooFormulaVar fL("fL","1 - fB - fC", RooArgList(fB,fC));



  RooAddPdf model("model","model",RooArgList(pdfB,pdfC,pdfL),RooArgList(fB,fC));

  //RooProdPdf modelConstrained("modelConstrained","",RooArgSet(model,cConstraint));

  // --------------------------
  // Fit
  // --------------------------

  // RooFitResult* res = model.fitTo(dataR,
  //                                 Save(),
  //                                 SumW2Error(kTRUE),
  //                                 PrintLevel(1));

  // res->Print();

  RooFitResult* res = model.fitTo(dataR,Save(),SumW2Error(kTRUE),PrintLevel(1));

  res->Print();

    res->correlationMatrix().Print();
  // --------------------------
  // Print fractions
  // --------------------------

  double fracb = fB.getVal();
  double fracc = fC.getVal();
  double fracl = 1.0 - fracb - fracc;

  printf("\nFractions:\n");
  printf("B = %.4f\n", fracb);
  printf("C = %.4f\n", fracc);
  printf("L = %.4f\n", fracl);

  Printf("Abundance Ratios : %f : %f : %f\n", fracb/fracb, fracc/fracb, fracl/fracb);

  TH1D *fractions = new TH1D("fractions", "", 3, 0, 3);
  fractions->SetBinContent(1, fL.getVal());
  fractions->SetBinError(1, TMath::Sqrt(TMath::Power(fB.getError(),2.) + TMath::Power(fC.getError(),2.)));
  
  fractions->SetBinContent(2, fC.getVal());
  fractions->SetBinError(2, fC.getError());

  fractions->SetBinContent(3, fB.getVal());
  fractions->SetBinError(3, fB.getError());

  TH2* hModel = (TH2*)hData->Clone("hModel");
hModel->Reset();

double nData = hData->Integral();

for(int ix=1; ix<=hData->GetNbinsX(); ++ix){
  for(int iy=1; iy<=hData->GetNbinsY(); ++iy){

    double xcenter = hData->GetXaxis()->GetBinCenter(ix);
    double ycenter = hData->GetYaxis()->GetBinCenter(iy);

    x.setVal(xcenter);
    y.setVal(ycenter);

    double val = model.getVal(obsSet);

    hModel->SetBinContent(ix, iy, val);
  }
}

hModel->Scale(1./hModel->Integral());

    TH1 *hxModel = (TH1*)hModel->ProjectionX("hxModel");
    TH1 *hxRatio = (TH1*)hxModel->Clone("hxRatio");
    TH1 *hxData = (TH1*)hData->ProjectionX("hxData");
    hxRatio->Divide(hxData, hxModel, 1., 1., "B");

    TH1 *hyModel = (TH1*)hModel->ProjectionY("hyModel");
    TH1 *hyRatio = (TH1*)hyModel->Clone("hyRatio");
    TH1 *hyData = (TH1*)hData->ProjectionY("hyData");
    hyRatio->Divide(hyData, hyModel, 1., 1., "B");
//gStyle->SetOptStat(0);
  // --------------------------
// X projection
// --------------------------

    RooPlot *xframe = x.frame();
    dataR.plotOn(xframe, Name("Data_x"));
    model.plotOn(xframe, LineColor(kRed), LineWidth(2), Name("model_fit_x"));
    model.plotOn(xframe, Components("pdfB"), LineColor(kBlue), LineStyle(kDashed), Name("bclf_b_plot"));
    model.plotOn(xframe, Components("pdfC"), LineColor(kGreen+2), LineStyle(kDashed), Name("bclf_c_plot"));
    model.plotOn(xframe, Components("pdfL"), LineColor(kMagenta), LineStyle(kDashed), Name("bclf_lf_plot"));

    RooHist *hist_data_x = (RooHist *) xframe->findObject("Data_x");
    RooCurve *bclf_c_fit = (RooCurve*) xframe->findObject("model_fit_x");
    RooCurve *bclf_c_b = (RooCurve*) xframe->findObject("bclf_b_plot");
    RooCurve *bclf_c_c = (RooCurve*) xframe->findObject("bclf_c_plot");
    RooCurve *bclf_c_lf = (RooCurve*) xframe->findObject("bclf_lf_plot");

    TLine *line = new TLine(0.,1., 1., 1.);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->SetLineStyle(7);

  TCanvas *cx = new TCanvas("cx","template fit",800,900);
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
  xframe->GetYaxis()->SetRangeUser(7e-6, 2e-1);
  xframe->SetTitle("");
  xframe->GetYaxis()->SetLabelSize(0.05);
  xframe->GetYaxis()->SetTitleSize(0.05);
  xframe->GetXaxis()->SetTitle("BDT(bc|usdg) score");
  xframe->GetYaxis()->SetTitle("Probability Density");

xframe->Draw();

TLegend *legendx = new TLegend(0.3,0.1,0.58,0.3);

    legendx->AddEntry(bclf_c_b, "b jets", "L");
    legendx->AddEntry(bclf_c_c, "c jets", "L");
    legendx->AddEntry(bclf_c_lf, "lf jets", "L");
    legendx->AddEntry(hist_data_x, "Data", "LEP");
    legendx->AddEntry(bclf_c_fit, "Fit", "L");
    legendx->SetBorderSize(0);
    legendx->SetFillColor(0);

legendx->Draw();

TPaveText *expt = new TPaveText(0.5,0.7,0.95,0.89,"NDC");
expt->AddText("pp, #sqrt{#it{s}} = 13.6 TeV");
expt->AddText("Anti-#it{k}_{T}, R = 0.4, |#eta| < 0.5");
expt->AddText(Form("%s GeV/#it{c}",ptbin.Data()));

expt->SetBorderSize(0);
expt->SetFillStyle(0);
expt->SetTextFont(42);

expt->Draw();

cx->Update();
pad2->cd();
  hxRatio->GetYaxis()->SetRangeUser(0.55, 1.45);
  hxRatio->GetYaxis()->SetLabelSize(0.1);
  hxRatio->GetYaxis()->SetTitleSize(0.12);
  hxRatio->GetYaxis()->SetTitleOffset(0.6);
  hxRatio->GetXaxis()->SetLabelSize(0.12);
  hxRatio->GetXaxis()->SetTitleSize(0.15);
  hxRatio->GetXaxis()->SetTickSize(0.05);
  hxRatio->GetYaxis()->SetNdivisions(7,5,0);
  hxRatio->GetXaxis()->SetTitle("BDT(bc|usdg) score");
  hxRatio->GetYaxis()->SetTitle("Data/Fit");
  hxRatio->SetMarkerStyle(kFullCircle);
  hxRatio->SetMarkerColor(kBlack);
  hxRatio->SetLineColor(kBlack);
  hxRatio->Draw();
  line->Draw();

cx->SaveAs("projection_x.pdf");


    RooPlot *yframe = y.frame();
    dataR.plotOn(yframe, Name("Data_y"));
    model.plotOn(yframe, LineColor(kRed), LineWidth(2), Name("model_fit_y"));
    model.plotOn(yframe, Components("pdfB"), LineColor(kBlue), LineStyle(kDashed), Name("bc_b_plot"));
    model.plotOn(yframe, Components("pdfC"), LineColor(kGreen+2), LineStyle(kDashed), Name("bc_c_plot"));
    model.plotOn(yframe, Components("pdfL"), LineColor(kMagenta), LineStyle(kDashed), Name("bc_lf_plot"));

    RooHist *hist_data_y = (RooHist *) yframe->findObject("Data_y");
    RooCurve *bc_c_fit = (RooCurve*) yframe->findObject("model_fit_y");
    RooCurve *bc_c_b = (RooCurve*) yframe->findObject("bc_b_plot");
    RooCurve *bc_c_c = (RooCurve*) yframe->findObject("bc_c_plot");
    RooCurve *bc_c_lf = (RooCurve*) yframe->findObject("bc_lf_plot");

  TCanvas *cy = new TCanvas("cy","template fit",800,900);
    gStyle->SetOptStat(0000);
  TPad *pad3 = new TPad("pad3","Main Pad",0,1-topPadFrac,1,1);
  TPad *pad4 = new TPad("pad4","Ratio Pad",0,0,1,1-topPadFrac);
  cy->cd();
  pad3->Draw();
  pad4->Draw();


  pad3->SetLeftMargin(leftMargin);
  pad3->SetRightMargin(rightMargin);
  pad3->SetTopMargin(topMargin);
  pad3->SetBottomMargin(0.001);
  pad3->SetLogy();
  pad3->SetTicks(1);

  pad4->SetLeftMargin(leftMargin);
  pad4->SetRightMargin(rightMargin);
  pad4->SetTopMargin(0.001);
  pad4->SetBottomMargin(bottomMargin);
  pad4->SetTicks(1);

  pad3->cd();
  yframe->GetYaxis()->SetRangeUser(7e-6, 2e-1);
  yframe->SetTitle("");
  yframe->GetYaxis()->SetLabelSize(0.05);
  yframe->GetYaxis()->SetTitleSize(0.05);
  yframe->GetXaxis()->SetTitle("BDT(b|c) score");
  yframe->GetYaxis()->SetTitle("Probability Density");

yframe->Draw();

TLegend *legendy = new TLegend(0.3,0.1,0.58,0.3);


    legendy->AddEntry(bc_c_b, "b jets", "L");
    legendy->AddEntry(bc_c_c, "c jets", "L");
    legendy->AddEntry(bc_c_lf, "lf jets", "L");
    legendy->AddEntry(hist_data_y, "Data", "LEP");
    legendy->AddEntry(bc_c_fit, "Fit", "L");
    legendy->SetBorderSize(0);
    legendy->SetFillColor(0);

legendy->Draw();

expt->Draw();

cy->Update();

pad4->cd();
  hyRatio->GetYaxis()->SetRangeUser(0.55, 1.45);
  hyRatio->GetYaxis()->SetLabelSize(0.1);
  hyRatio->GetYaxis()->SetTitleSize(0.12);
  hyRatio->GetYaxis()->SetTitleOffset(0.6);
  hyRatio->GetXaxis()->SetLabelSize(0.12);
  hyRatio->GetXaxis()->SetTitleSize(0.15);
  hyRatio->GetXaxis()->SetTickSize(0.05);
  hyRatio->GetYaxis()->SetNdivisions(7,5,0);
  hyRatio->GetXaxis()->SetTitle("BDT(b|c) score");
  hyRatio->GetYaxis()->SetTitle("Data/Fit");
  hyRatio->SetMarkerStyle(kFullCircle);
  hyRatio->SetMarkerColor(kBlack);
  hyRatio->SetLineColor(kBlack);
  hyRatio->Draw();
  line->Draw();
cy->SaveAs("projection_y.pdf");


// Draw comparison
TCanvas* c3 = new TCanvas("c3","2D comparison",800,900);
c3->Divide(2,2);

c3->cd(1);
gPad->SetLogz();
gPad->SetRightMargin(0.15);
hData->SetTitle("Data");
hData->GetXaxis()->SetTitle("BDT(bc|lf)");
hData->GetYaxis()->SetTitle("BDT(b|c)");
hData->Draw("COLZ");
gPad->SetTitle("Data");

c3->cd(2);
gPad->SetLogz();
gPad->SetRightMargin(0.15);
hModel->SetMinimum(0);
hModel->SetTitle("Model");
hModel->GetXaxis()->SetTitle("BDT(bc|lf)");
hModel->GetYaxis()->SetTitle("BDT(b|c)");
hModel->Draw("COLZ");
gPad->SetTitle("Model");

c3->cd(3);

gPad->SetLogz();
gPad->SetRightMargin(0.15);
TH2* hResidual = (TH2*)hData->Clone("hResidual");
hResidual->Add(hModel, -1);
hResidual->GetXaxis()->SetTitle("BDT(bc|lf)");
hResidual->GetYaxis()->SetTitle("BDT(b|c)");
hResidual->SetTitle("Data - Model");
hResidual->Draw("COLZ");

c3->cd(4);
TH2* hPullMap = (TH2*)hData->Clone("hPullMap");

for(int ix=1; ix<=hData->GetNbinsX(); ++ix){
  for(int iy=1; iy<=hData->GetNbinsY(); ++iy){

    double d = hData->GetBinContent(ix,iy);
    double m = hModel->GetBinContent(ix,iy);
    double sd = hData->GetBinError(ix,iy);
    double sm = hModel->GetBinError(ix,iy);
    double s = TMath::Sqrt(TMath::Power(sd,2.) + TMath::Power(sm,2.));

    if(m > 0){
      hPullMap->SetBinContent(ix,iy, (d-m)/s);
    } else {
      hPullMap->SetBinContent(ix,iy, 0);
    }
  }
}
//gPad->SetLogz();
hPullMap->SetTitle("Pull Map");
hPullMap->Draw("Colz");

c3->SaveAs("2d_comparison.pdf");

//Plot the NEGATIVE LOG LIKELIHOOD

// TCanvas* cnll = new TCanvas("cnll","NLL surface",1500,400);
// cnll->Divide(3,1);
// cnll->cd(1);
// RooAbsReal* nllB = model.createNLL(dataR);
// RooAbsReal* pllB = nllB->createProfile(RooArgSet(fB));
// RooPlot* frameB = fB.frame();
// pllB->plotOn(frameB,ShiftToZero());
// frameB->Draw();

// cnll->cd(2);
// RooAbsReal* nllC = model.createNLL(dataR);
// RooAbsReal* pllC = nllC->createProfile(RooArgSet(fC));
// RooPlot* frameC = fC.frame();
// pllC->plotOn(frameC,ShiftToZero());
// frameC->Draw();

// cnll->cd(3);
// RooAbsReal* nllL = model.createNLL(dataR);
// RooAbsReal* pllL = nllL->createProfile(RooArgSet(fL));
// RooPlot* frameL = fL.frame();
// pllL->plotOn(frameL,ShiftToZero());
// frameL->Draw();

  TFile *fitresult = new TFile("fitresult2D.root", "recreate");
  fitresult->cd();

  res->Write();
  cx->Write();
  cy->Write();
  c3->Write();
  hData->Write();
  hModel->Write();
  hResidual->Write();
  fractions->Write();
  hxRatio->Write();
  hyRatio->Write();
  fitresult->Write();
  fitresult->Close();
  

}