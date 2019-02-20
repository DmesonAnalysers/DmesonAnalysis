//___________________________________________________________________________________//
// Macro for fitting D+ and Ds+ invariant-mass spectra                               //
// Main Function: GetRawYieldsDplusDs                                                //
// Author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino                    //
//___________________________________________________________________________________//

#if !defined (__CINT__) || defined (__CLING__)

#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include "TROOT.h"
#include "Riostream.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TDatabasePDG.h"

#include "AliHFInvMassFitter.h"
#include "AliVertexingHFUtils.h"

#endif

using namespace std;

enum {k010, k3050};
enum {kDplus, kDs};

//__________________________________________________________________________________________________________________
int GetRawYieldsDplusDs(int cent = k010, bool isMC = kFALSE, TString infilename = "InvMassSpectraDplus_010_PbPb2015cuts.root", TString cfgfilename = "Dplus/config_Dplus_Fit.yml", TString outFileName = "RawYieldsDplus_010_PbPb2015cuts.root");
double SingleGaus(double *m, double *pars);
double DoublePeakSingleGaus(double *x, double *pars);
void SetHistoStyle(TH1 *histo, int color=kBlack);
void SetStyle();
void DivideCanvas(TCanvas* c, int nPtBins);

//__________________________________________________________________________________________________________________
int GetRawYieldsDplusDs(int cent, bool isMC, TString infilename, TString cfgfilename, TString outFileName) {

  SetStyle();

  //load config
  TString centname = "";
  if(cent==k010) centname = "Cent010";
  else if(cent==k3050) centname = "Cent3050";

  YAML::Node config = YAML::LoadFile(cfgfilename.Data());
  string mesonname = config[centname.Data()]["Meson"].as<string>();
  int meson = (mesonname=="Dplus") ? kDplus : kDs;
  bool fixSigma = config[centname.Data()]["FixSigma"].as<int>();
  bool fixMean = config[centname.Data()]["FixMean"].as<int>();
  bool UseLikelihood = config[centname.Data()]["UseLikelihood"].as<int>();
  vector<double> PtMin = config[centname.Data()]["PtMin"].as<vector<double>>();
  vector<double> PtMax = config[centname.Data()]["PtMax"].as<vector<double>>();
  vector<double> MassMin = config[centname.Data()]["MassMin"].as<vector<double>>();
  vector<double> MassMax = config[centname.Data()]["MassMax"].as<vector<double>>();
  vector<int> Rebin = config[centname.Data()]["Rebin"].as<vector<int>>();
  vector<int> InclSecPeak = config[centname.Data()]["InclSecPeak"].as<vector<int>>();
  vector<string> bkgfunc = config[centname.Data()]["BkgFunc"].as<vector<string>>();
  vector<string> sgnfunc = config[centname.Data()]["SgnFunc"].as<vector<string>>();
  const unsigned int nPtBins = PtMin.size();
  int BkgFunc[nPtBins], SgnFunc[nPtBins];
  double PtLims[nPtBins+1];

  for(int iPt=0; iPt<nPtBins; iPt++) {
    PtLims[iPt] = PtMin[iPt];
    PtLims[iPt+1] = PtMax[iPt];
    if(bkgfunc[iPt] == "kExpo") BkgFunc[iPt] = AliHFInvMassFitter::kExpo; 
    else if(bkgfunc[iPt] == "kLin") BkgFunc[iPt] = AliHFInvMassFitter::kLin; 
    else if(bkgfunc[iPt] == "kPol2") BkgFunc[iPt] = AliHFInvMassFitter::kPol2; 
    if(sgnfunc[iPt] == "kGaus") SgnFunc[iPt] = AliHFInvMassFitter::kGaus; 
    else if(sgnfunc[iPt] == "k2Gaus") SgnFunc[iPt] = AliHFInvMassFitter::k2Gaus; 
  }

  //load inv-mass histos
  auto infile = TFile::Open(infilename.Data());
  if(!infile || !infile->IsOpen()) return -1;
  TH1F* hMass[nPtBins];
  TH1F* hEv = NULL;
  for(int iPt=0; iPt<nPtBins; iPt++) {
    hMass[iPt] = static_cast<TH1F*>(infile->Get(Form("hMass_%0.f_%0.f",PtMin[iPt],PtMax[iPt])));
    hEv = static_cast<TH1F*>(infile->Get("hEvForNorm"));
    hMass[iPt]->SetDirectory(0);
    hEv->SetDirectory(0);
    SetHistoStyle(hMass[iPt]);
    SetHistoStyle(hEv);
  }
  infile->Close();

  //define output histos
  auto hRawYields = new TH1D("hRawYields",";#it{p}_{T} (GeV/#it{c});raw yield",nPtBins,PtLims);
  auto hRawYieldsSigma = new TH1D("hRawYieldsSigma",";#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})",nPtBins,PtLims);
  auto hRawYieldsSigma2 = new TH1D("hRawYieldsSigma2",";#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})",nPtBins,PtLims);
  auto hRawYieldsMean = new TH1D("hRawYieldsMean",";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})",nPtBins,PtLims);
  auto hRawYieldsFracGaus2 = new TH1D("hRawYieldsFracGaus2",";#it{p}_{T} (GeV/#it{c});second-gaussian fraction",nPtBins,PtLims);
  auto hRawYieldsSignificance = new TH1D("hRawYieldsSignificance",";#it{p}_{T} (GeV/#it{c});significance (3#sigma)",nPtBins,PtLims);
  auto hRawYieldsSoverB = new TH1D("hRawYieldsSoverB",";#it{p}_{T} (GeV/#it{c});S/B (3#sigma)",nPtBins,PtLims);
  auto hRawYieldsSignal = new TH1D("hRawYieldsSignal",";#it{p}_{T} (GeV/#it{c});Signal (3#sigma)",nPtBins,PtLims);
  auto hRawYieldsBkg = new TH1D("hRawYieldsBkg",";#it{p}_{T} (GeV/#it{c});Background (3#sigma)",nPtBins,PtLims);
  auto hRawYieldsChiSquare = new TH1D("hRawYieldsChiSquare",";#it{p}_{T} (GeV/#it{c});#chi^{2}/#it{ndf}",nPtBins,PtLims);
  auto hRawYieldsSecondPeak = new TH1D("hRawYieldsSecondPeak",";#it{p}_{T} (GeV/#it{c});raw yield second peak",nPtBins,PtLims);
  auto hRawYieldsMeanSecondPeak = new TH1D("hRawYieldsMeanSecondPeak",";#it{p}_{T} (GeV/#it{c});mean second peak (GeV/#it{c}^{2})",nPtBins,PtLims);
  auto hRawYieldsSigmaSecondPeak = new TH1D("hRawYieldsSigmaSecondPeak",";#it{p}_{T} (GeV/#it{c});width second peak (GeV/#it{c}^{2})",nPtBins,PtLims);
  auto hRawYieldsSigma2SecondPeak = new TH1D("hRawYieldsSigma2SecondPeak",";#it{p}_{T} (GeV/#it{c});width second peak (GeV/#it{c}^{2})",nPtBins,PtLims);
  auto hRawYieldsFracGaus2SecondPeak = new TH1D("hRawYieldsFracGaus2SecondPeak",";#it{p}_{T} (GeV/#it{c});second-gaussian fraction second peak (GeV/#it{c}^{2})",nPtBins,PtLims);
  auto hRawYieldsSignificanceSecondPeak = new TH1D("hRawYieldsSignificanceSecondPeak",";#it{p}_{T} (GeV/#it{c});signficance second peak (3#sigma)",nPtBins,PtLims);
  auto hRawYieldsSigmaRatioSecondFirstPeak = new TH1D("hRawYieldsSigmaRatioSecondFirstPeak",";#it{p}_{T} (GeV/#it{c});width second peak / width first peak",nPtBins,PtLims);
  auto hRawYieldsTrue = new TH1D("hRawYieldsTrue",";#it{p}_{T} (GeV/#it{c});true signal",nPtBins,PtLims);
  auto hRawYieldsSecondPeatrue = new TH1D("hRawYieldsSecondPeatrue",";#it{p}_{T} (GeV/#it{c});true signal second peak",nPtBins,PtLims);
  auto hRelDiffRawYieldsFitTrue = new TH1D("hRelDiffRawYieldsFitTrue",";#it{p}_{T} (GeV/#it{c}); (Y_{fit} - Y_{true}) / Y_{true}",nPtBins,PtLims);
  auto hRelDiffRawYieldsSecondPeakFitTrue = new TH1D("hRelDiffRawYieldsSecondPeakFitTrue",";#it{p}_{T} (GeV/#it{c});(Y_{fit} - Y_{true}) / Y_{true} second peak",nPtBins,PtLims);
  SetHistoStyle(hRawYields);
  SetHistoStyle(hRawYieldsSigma);
  SetHistoStyle(hRawYieldsSigma2);
  SetHistoStyle(hRawYieldsMean);
  SetHistoStyle(hRawYieldsFracGaus2);
  SetHistoStyle(hRawYieldsSignificance);
  SetHistoStyle(hRawYieldsSoverB);
  SetHistoStyle(hRawYieldsSignal);
  SetHistoStyle(hRawYieldsBkg);
  SetHistoStyle(hRawYieldsChiSquare);
  SetHistoStyle(hRawYieldsSecondPeak,kRed+1);
  SetHistoStyle(hRawYieldsMeanSecondPeak,kRed+1);
  SetHistoStyle(hRawYieldsSigmaSecondPeak,kRed+1);
  SetHistoStyle(hRawYieldsSigma2SecondPeak,kRed+1);
  SetHistoStyle(hRawYieldsFracGaus2SecondPeak,kRed+1);
  SetHistoStyle(hRawYieldsSignificanceSecondPeak,kRed+1);
  SetHistoStyle(hRawYieldsSigmaRatioSecondFirstPeak,kRed+1);
  SetHistoStyle(hRawYieldsTrue);
  SetHistoStyle(hRawYieldsSecondPeatrue,kRed+1);
  SetHistoStyle(hRelDiffRawYieldsFitTrue);
  SetHistoStyle(hRelDiffRawYieldsSecondPeakFitTrue,kRed+1);
  
  //fit histos
  TCanvas* cMass = new TCanvas("cMass","Mass",1920,1080);
  DivideCanvas(cMass,nPtBins);

  double massDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  double massDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
  double massForFit = (meson==kDplus) ? massDplus : massDs;

  for(int iPt=0; iPt<nPtBins; iPt++) {

    auto hMassForFit=(TH1F*)AliVertexingHFUtils::RebinHisto(hMass[iPt],Rebin[iPt]);
    hMassForFit->SetName(Form("MassForFit%d",iPt));
    SetHistoStyle(hMassForFit);

    if(isMC) {
      TF1* massFunc = NULL;
      if(SgnFunc[iPt]==AliHFInvMassFitter::kGaus) {
        if(!InclSecPeak[iPt]) massFunc = new TF1(Form("massFunc%d",iPt),SingleGaus,MassMin[iPt],MassMax[iPt],3); 
        else massFunc = new TF1(Form("massFunc%d",iPt),DoublePeakSingleGaus,MassMin[iPt],MassMax[iPt],6); 
      }
      else {
        cout << "To be implemented" << endl;
        return -1;
      }
    }
    else {
      auto massFitter = new AliHFInvMassFitter(hMassForFit,MassMin[iPt],MassMax[iPt],BkgFunc[iPt],SgnFunc[iPt]);
      if(UseLikelihood) massFitter->SetUseLikelihoodFit();
      massFitter->SetInitialGaussianMean(massForFit);
      massFitter->SetInitialGaussianSigma(0.010);
      if(InclSecPeak[iPt] && meson==kDs) massFitter->IncludeSecondGausPeak(massDplus,false,0.008,true);
      massFitter->MassFitter(false);
      
      double rawyield = massFitter->GetRawYield();
      double rawyielderr = massFitter->GetRawYieldError();
      double sigma = massFitter->GetSigma();
      double sigmaerr = massFitter->GetSigmaUncertainty();
      double mean = massFitter->GetMean();
      double meanerr = massFitter->GetMeanUncertainty();
      double redchi2 = massFitter->GetReducedChiSquare();
      double signif=0., signiferr=0.;
      double sgn=0., sgnerr=0.;
      double bkg=0., bkgerr=0.;
      massFitter->Significance(3,signif,signiferr);
      massFitter->Signal(3,sgn,sgnerr);
      massFitter->Background(3,bkg,bkgerr);

      hRawYields->SetBinContent(iPt+1,rawyield);
      hRawYields->SetBinError(iPt+1,rawyielderr);
      hRawYieldsSigma->SetBinContent(iPt+1,sigma);
      hRawYieldsSigma->SetBinError(iPt+1,sigmaerr);
      hRawYieldsMean->SetBinContent(iPt+1,mean);
      hRawYieldsMean->SetBinError(iPt+1,meanerr);
      hRawYieldsSignificance->SetBinContent(iPt+1,signif);
      hRawYieldsSignificance->SetBinError(iPt+1,signiferr);
      hRawYieldsSoverB->SetBinContent(iPt+1,sgn/bkg);
      hRawYieldsSoverB->SetBinError(iPt+1,sgn/bkg*TMath::Sqrt(sgnerr/sgn*sgnerr/sgn+bkgerr/bkg*bkgerr/bkg));
      hRawYieldsSignal->SetBinContent(iPt+1,sgn);
      hRawYieldsSignal->SetBinError(iPt+1,sgnerr);
      hRawYieldsBkg->SetBinContent(iPt+1,bkg);
      hRawYieldsBkg->SetBinError(iPt+1,bkgerr);
      hRawYieldsChiSquare->SetBinContent(iPt+1,redchi2);
      hRawYieldsChiSquare->SetBinError(iPt+1,1.e-20);
      
      if(nPtBins>1) 
        cMass->cd(iPt+1);
      hMassForFit->GetYaxis()->SetRangeUser(hMassForFit->GetMinimum()*0.95,hMassForFit->GetMaximum()*1.2);
      massFitter->DrawHere(gPad);
    }
    cMass->Modified();
    cMass->Update();
  }

  //save output histos
  TFile outFile(outFileName.Data(),"RECREATE");
  cMass->Write();
  for(int iPt=0; iPt<nPtBins; iPt++)
    hMass[iPt]->Write();
  hRawYields->Write();
  hRawYieldsSigma->Write();
  hRawYieldsMean->Write();
  hRawYieldsSignificance->Write();
  hRawYieldsSoverB->Write();
  hRawYieldsSignal->Write();
  hRawYieldsBkg->Write();
  hRawYieldsChiSquare->Write();
  hEv->Write();
  outFile.Close();

  outFileName.ReplaceAll(".root",".pdf");
  cMass->SaveAs(outFileName.Data());

  return 0;
}

//__________________________________________________________________________________________________________________
double SingleGaus(double *m, double *pars) {
  double norm = pars[0], mean = pars[1], sigma = pars[2];

  return norm*TMath::Gaus(m[0],mean,sigma,true);
}

//__________________________________________________________________________________________________________________
double DoublePeakSingleGaus(double *x, double *pars) {
  double norm1 = pars[0], mean1 = pars[1], sigma1 = pars[2];
  double norm2 = pars[3], mean2 = pars[4], sigma2 = pars[5];

  return norm1*TMath::Gaus(x[0],mean1,sigma1,true) + norm2*TMath::Gaus(x[0],mean2,sigma2,true);
}

//__________________________________________________________________________________________________________________
void SetHistoStyle(TH1 *histo, int color) {
  histo->SetStats(kFALSE);
  histo->SetMarkerSize(1.);
  histo->SetMarkerStyle(20);
  histo->SetLineWidth(2);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}

//__________________________________________________________________________________________________________________
void DivideCanvas(TCanvas* c, int nPtBins) {
  if(nPtBins<2)
    c->cd();
  else if(nPtBins==2 || nPtBins==3)
    c->Divide(nPtBins,1);
  else if(nPtBins==4 || nPtBins==6 || nPtBins==8)
    c->Divide(nPtBins/2,2);
  else if(nPtBins==5 || nPtBins==7)
    c->Divide((nPtBins+1)/2,2);
  else if(nPtBins==9 || nPtBins==12 || nPtBins==15)
    c->Divide(nPtBins/3,3);
  else if(nPtBins==10 || nPtBins==11)
    c->Divide(4,3);
  else if(nPtBins==13 || nPtBins==14)
    c->Divide(5,3);
  else if(nPtBins>15 && nPtBins<=20 && nPtBins%4==0)
    c->Divide(nPtBins/4,4);
  else if(nPtBins>15 && nPtBins<=20 && nPtBins%4!=0)
    c->Divide(5,4);
  else if(nPtBins==21)
    c->Divide(7,3);
  else if(nPtBins>21 && nPtBins<=25)
    c->Divide(5,5);
  else if(nPtBins>25 && nPtBins%2==0)
    c->Divide(nPtBins/2,2);
  else
    c->Divide((nPtBins+1)/2,2);
}

//__________________________________________________________________________________________________________________
void SetStyle() {
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TGaxis::SetMaxDigits(3);
}
