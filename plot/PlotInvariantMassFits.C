#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TASImage.h>
#include <TPad.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <Riostream.h>
#include <TColor.h>
#include <TGaxis.h>
#endif

const int markerstyle=kFullCircle;
const double markersize=.6;
const int markercolor=kBlack;
const int histolinecolor=kBlack;
const int totfunccolor=kBlue;
const int bkgfunccolor=kRed;
const int totfuncstyle=1;
const int bkgfuncstyle=2;
const double funcwidth=2;
const int npx=200;

const double bottommargin=0.12;
const double leftmargin=0.15;
const double rightmargin=0.015;
const double topmargin=0.06;
const double xaxisoffset=1.;
const double yaxisoffset=1.1;
const int maxdigits=3;

const int titlefont=42;
const double titlesize=0.04;
const int labelfont=42;
const double labelsize=0.03;

const int textfont=42;
const double textsize=0.04;
const double textsizeLab=0.04;

double lowEdgeOfThirdNonZeroBin(TH1F *histo){
  for(int i = 1; i <= histo->GetNbinsX(); i++){
    if(histo->GetBinContent(i) > 0){
      return histo->GetBinLowEdge(i+2);
    }
  }
  printf("Warning: all bins are zero");
  return nan("");
}

double lowEdgeOfPenultimateNonZeroBin(TH1F *histo){
  for(int i = histo->GetNbinsX(); i >= 1; i--){
    if(histo->GetBinContent(i) > 0){
      return histo->GetBinLowEdge(i-1);
    }
  }
  printf("Warning: all bins are zero");
  return nan("");
}

int PlotInvariantMassFitsLc() {
  TString residuals;
  residuals = "_Residuals";
  residuals = "";

  TString specifier="FDEn";
  TString infilename=Form("rawYields/basic/ry_basic_%s.root", specifier.Data());
  TString canvasname="cResiduals0";
  TString histoname="fHistoInvMass";
  TString totfuncname="funcmass";
  TString bkgfuncname="funcbkgrefit";

  TGaxis::SetMaxDigits(maxdigits);
  gStyle->SetPadBottomMargin(bottommargin);
  gStyle->SetPadLeftMargin(leftmargin);
  gStyle->SetPadRightMargin(rightmargin);
  gStyle->SetPadTopMargin(topmargin);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetTitleSize(textsize,"xyz");
  gStyle->SetLabelSize(textsize,"xyz");
  //gStyle->SetNdivisions(505);

  TString particlename="#Lambda_{c}^{+} #rightarrow pK^{0}_{s}";
  TString massaxistit="#it{M}(pK^{0}_{s}) (GeV/#it{c}^{2})";

  TH1F* hMass[4];
  TF1* fTot[4];
  TF1* fBkg[4];

  TString widthtext[4];
  TString meantext[4];
  TString signaltext[4];
  TString signiftext[4];
  TString SoverBtext[4];

  TString pttitle[4] = {"2 < #it{p}_{T} < 4 GeV/#it{c}",
                        "4 < #it{p}_{T} < 6 GeV/#it{c}",
                        "6 < #it{p}_{T} < 8 GeV/#it{c}",
                        "8 < #it{p}_{T} < 12 GeV/#it{c}"};
  int ptbin[4] = {1, 2, 3, 4};
  TFile* infile = TFile::Open(infilename.Data(),"READ");
  if(!infile) {return 1;}
  TCanvas* cinput=(TCanvas*)infile->Get(canvasname.Data());
  if(!cinput) {cerr << "Error: wrong canvas name, check it!" <<endl; return 2;}
  for(int iPt=0; iPt<4; iPt++) {
    hMass[iPt] = (TH1F*)cinput->GetPad(ptbin[iPt])->GetPrimitive(histoname.Data());
    if(!hMass[iPt]) {cerr << "Error: wrong histo name, check it!" <<endl; return iPt+1;}
    hMass[iPt]->SetMarkerStyle(markerstyle);
    hMass[iPt]->SetMarkerSize(markersize);
    hMass[iPt]->SetMarkerColor(markercolor);
    hMass[iPt]->SetLineColor(histolinecolor);
    hMass[iPt]->GetXaxis()->SetTitle(massaxistit.Data());
    hMass[iPt]->GetYaxis()->SetTitle(Form("Counts per %0.f MeV/#it{c}^{2}",hMass[iPt]->GetBinWidth(5)*1000));
    hMass[iPt]->GetXaxis()->SetTitleOffset(xaxisoffset);
    hMass[iPt]->GetYaxis()->SetTitleOffset(yaxisoffset);
    hMass[iPt]->GetXaxis()->SetTitleSize(titlesize);
    hMass[iPt]->GetYaxis()->SetTitleSize(titlesize);
    hMass[iPt]->GetXaxis()->SetTitleFont(titlefont);
    hMass[iPt]->GetYaxis()->SetTitleFont(titlefont);
    hMass[iPt]->GetXaxis()->SetLabelSize(labelsize);
    hMass[iPt]->GetYaxis()->SetLabelSize(labelsize);
    hMass[iPt]->GetXaxis()->SetLabelFont(labelfont);
    hMass[iPt]->GetYaxis()->SetLabelFont(labelfont);
    hMass[iPt]->GetYaxis()->SetDecimals();
    hMass[iPt]->SetTitle("");
    fTot[iPt] = (TF1*)cinput->GetPad(ptbin[iPt])->GetPrimitive(totfuncname.Data());
    fTot[iPt]->SetLineColor(totfunccolor);
    fTot[iPt]->SetLineWidth(funcwidth);
    fTot[iPt]->SetLineStyle(totfuncstyle);
    fTot[iPt]->SetNpx(npx);
    fBkg[iPt] = (TF1*)cinput->GetPad(ptbin[iPt])->GetPrimitive(bkgfuncname.Data());
    fBkg[iPt]->SetLineColor(bkgfunccolor);
    fBkg[iPt]->SetLineWidth(funcwidth);
    fBkg[iPt]->SetLineStyle(bkgfuncstyle);
    fBkg[iPt]->SetNpx(npx);
  

    double mean = fTot[iPt]->GetParameter(4);
    double sigma = fTot[iPt]->GetParameter(5);
    double errmean = fTot[iPt]->GetParError(4);
    double errsigma = fTot[iPt]->GetParError(5);
    double rawyield = fTot[iPt]->GetParameter(3)/hMass[iPt]->GetBinWidth(1);
    double rawyielderr = fTot[iPt]->GetParError(3)/hMass[iPt]->GetBinWidth(1);
    std::cout<<"mean: "<<mean<<"   sigma: "<<sigma<<"   rawyield: "<<rawyield<< std::endl;
    double ints=fTot[iPt]->Integral(mean-3*sigma,mean+3*sigma, 1e-5)/hMass[iPt]->GetBinWidth(4);
    double intb=fBkg[iPt]->Integral(mean-3*sigma,mean+3*sigma, 1e-5)/hMass[iPt]->GetBinWidth(2);
    double signal = ints-intb;
    double signalerr = rawyield/rawyielderr*signal;
    double bkg = intb;
    double bkgerr = fBkg[iPt]->GetParError(0)/fBkg[iPt]->GetParameter(0)*bkg;
    double significance = signal/TMath::Sqrt(signal+bkg);
    double significanceerr = significance*TMath::Sqrt((signalerr*signalerr+bkgerr*bkgerr)/(4.*(signal+bkg)*(signal+bkg))+(bkg/(signal+bkg))*(signalerr*signalerr)/signal/signal);
        double signaloverbkg = signal/bkg;

    if(iPt<10) {
      meantext[iPt]=Form("#mu = (%0.1f #pm %0.1f) MeV/#it{c}^{2}",mean*1000,errmean*1000);
      widthtext[iPt]=Form("#sigma = (%0.1f #pm %0.1f) MeV/#it{c}^{2}",sigma*1000,errsigma*1000);
    }
    else {
      meantext[iPt]=Form("#mu = (%0.1f #pm %0.1f) MeV/#it{c}^{2}",mean*1000,errmean*1000);
      widthtext[iPt]=Form("#sigma = %0.1f MeV/#it{c}^{2} fix to MB",sigma*1000);
    }
    signaltext[iPt]=Form("S = %0.f #pm %0.f",rawyield,rawyielderr);
    SoverBtext[iPt]=Form("S/B(3#sigma) = %0.3f",signaloverbkg);
  }


  TLatex* lat = new TLatex();
  lat->SetTextFont(textfont);
  lat->SetTextSize(textsizeLab);
  lat->SetNDC();

  TLatex* lat2 = new TLatex();
  lat2->SetTextFont(textfont);
  lat2->SetTextSize(textsize);
  lat2->SetNDC();

  TLegend* leg = new TLegend(0.22,0.7,0.58,0.9);
  leg->SetBorderSize(0);
  leg->SetTextSize(textsize);
  leg->SetFillStyle(0);
  leg->AddEntry(hMass[0],"Data","p");
  leg->AddEntry(fTot[0],"Total fit function","l");
  leg->AddEntry(fBkg[0],"Combinatorial background","l");

  // prompt mass
  // double ymin[4] = {120e3, 25e3, 2e3, 300};
  // double ymax[4] = {160e3, 31e3, 3.8e3, 950};
  // FD mass
  // double ymin[4] = {2800, 200, 0, 30};
  // double ymax[4] = {4500, 700, 120, 280};
  
  // prompt mass residuals
  // double ymin[4] = {120e3, 25e3, 2e3, 300};
  // double ymax[4] = {160e3, 31e3, 3.8e3, 950};
  // FD mass
  double ymin[4] = {0, 0, 0, 0};
  double ymax[4] = {1600, 300, 160, 100};


  double xinfo[4] = {0.22,0.22,0.22,0.22};
  double yinfo[4] = {0.35,0.35,0.35,0.35};
  double centmin[4] = {0,0,0,0};
  double centmax[4] = {10,10,10,10};
  double ytitle[4] = {0.44,0.44,0.44,0.44};
  int draworder[4] = {1,2,3,4};


  TCanvas *coutput = new TCanvas("coutput","",900,900);
  coutput->Divide(2,2);
  TH1F* hFrame[4];
  for(int iPt=0; iPt<4; iPt++) {
    double leftx = lowEdgeOfThirdNonZeroBin(hMass[iPt]);
    double rightx = lowEdgeOfPenultimateNonZeroBin(hMass[iPt]);

    hFrame[iPt] = coutput->cd(draworder[iPt])->DrawFrame(leftx,ymin[iPt],rightx,ymax[iPt],Form(";%s;%s",hMass[iPt]->GetXaxis()->GetTitle(),hMass[iPt]->GetYaxis()->GetTitle()));
    hFrame[iPt]->GetYaxis()->SetDecimals(2);

    TPaveText *pt = new TPaveText(leftmargin+0.1, 0.95, 1-rightmargin-0.1, 1, "NBNDC");
    pt->SetFillColor(0);
    pt->SetTextFont(textfont);
    pt->SetTextSize(textsize);
    pt->AddText(pttitle[iPt]);
    pt->Draw("");

    hMass[iPt]->Draw("E  same");
    fBkg[iPt]->Draw("same");
    fTot[iPt]->Draw("same");
    hMass[iPt]->GetYaxis()->SetRangeUser(ymin[iPt], ymax[iPt]);
    
  
    if(iPt==0) {
      //lat->DrawLatex(0.73,0.85,"This Thesis");
      lat2->DrawLatex(0.7,0.85,"pp, #sqrt{#it{s}} = 13 TeV");
      lat2->DrawLatex(0.7,0.75,particlename.Data());
    }
    else if(iPt==1) {
    }
    else if(iPt==0){
       leg->Draw("same");
    }
    //lat2->SetTextAlign(12);
    //lat2->DrawLatex((1 - rightmargin - leftmargin)/2, 0.97, pttitle[iPt].Data());
    

    lat2->DrawLatex(xinfo[iPt],yinfo[iPt],meantext[iPt].Data());
    lat2->DrawLatex(xinfo[iPt],yinfo[iPt]-0.06,widthtext[iPt].Data());
    lat2->DrawLatex(xinfo[iPt],yinfo[iPt]-0.12,signaltext[iPt].Data());
    lat2->DrawLatex(xinfo[iPt],yinfo[iPt]-0.18,SoverBtext[iPt].Data());
  }




  //output
  coutput->SaveAs(Form("MassLc_%s_%s.pdf", residuals.Data(), specifier.Data()));

  return 0;
}
