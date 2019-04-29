#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPad.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TGaxis.h>
#include <TPaveText.h>
#include <TArrayD.h>
#include <TNtuple.h>

#include "AliVertexingHFUtils.h"
#include "AliHFInvMassFitter.h"

#endif

//macro for the standard analysis of D+ mesons
//author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino

//*********************************************//
//                                             //
//    Main Function: RawYieldSystematics       //
//                                             //
//*********************************************//

//________________________________________________________________________________________________________________
//global variables
const TString reffilename = "outputs/rawyields/RawYieldsDs_010_central_2018.root";
const TString reffileMCname = "outputs/rawyields/RawYieldsDs_MC_010_central_2018.root";
const TString reffileMCname2 = "";
const int ptshifttoref = 0;

enum {kFree,kFixed,kFixedPlus15Perc,kFixedMinus15Perc};

const int nMins = 5;
const int nMaxs = 5;
const int nReb = 5;
const int nBkgFunc = 3;
const int nSgnFunc = 1;
const int nSigma = 2;
const int nMean = 1;
const double mins[nMins] = {1.78,1.79,1.80,1.81,1.82};
const double maxs[nMaxs] = {2.08,2.1,2.12,2.14,2.16};
const int rebin[nReb] = {2,3,4,5,6};
const int sgnfcn[nSgnFunc] = {AliHFInvMassFitter::kGaus};
const int bkgfcn[nBkgFunc] = {AliHFInvMassFitter::kExpo,AliHFInvMassFitter::kLin,AliHFInvMassFitter::kPol2};
const double maxchisquare = 2.;
const double minchisquare = 0.2;
const int nBinCounting = 2;
const double nSigmaBinCount[2] = {5,3};

//________________________________________________________________________________________________________________
//functions prototypes
int RawYieldSystematics(TString outfilerawname="outputs/rawyieldsyst/RawYieldsSystDs_010_central_2018.root");
double min(TH1F* histo);
double max(TH1F* histo);
int LoadRefFiles(TString reffilename, TString reffileMCname, TString reffileMCname2, TH1F *&hRawYieldRef, TH1F *&hSigmaRef, TH1F *&hMeanRef, TH1F *&hBkgRef, TH1F *&hSigmaMC, TH1F *&hFracGaus2MC, TH1F *&hSigmaGaus2MC);
void SetStyle();

//________________________________________________________________________________________________________________
int RawYieldSystematics(TString outfilerawname) {
    
  //setting drawing style
  SetStyle();
  
  //load histos
  TH1F* hRawYieldRef=0x0;
  TH1F* hSigmaRef=0x0;
  TH1F* hMeanRef=0x0;
  TH1F* hBkgRef=0x0;
  TH1F* hSigmaMC=0x0;
  TH1F* hFracGaus2MC=0x0;
  TH1F* hRatioSigmaGaus2MC=0x0;
  int loadref = LoadRefFiles(reffilename,reffileMCname,reffileMCname2,hRawYieldRef,hSigmaRef,hMeanRef,hBkgRef,hSigmaMC,hFracGaus2MC,hRatioSigmaGaus2MC);
  if(loadref>0 && loadref<9) {cerr << "Error in loading reference files! Check them please." <<endl; return loadref;}
  if(hRatioSigmaGaus2MC && hSigmaMC) {
    hRatioSigmaGaus2MC->Divide(hSigmaMC);
  }

  const int nPtBins = hRawYieldRef->GetNbinsX(); //2 values for each cut variable (min,max)
  const int nPtLims = nPtBins+1;
  double PtLims[nPtLims];
  TH1F* hMass[nPtBins];
  TFile* infile = TFile::Open(reffilename.Data());
  for(int iPt=0; iPt<nPtBins; iPt++) {
    PtLims[iPt] = hRawYieldRef->GetBinLowEdge(iPt+1);
  }
  PtLims[nPtBins] = hRawYieldRef->GetBinLowEdge(nPtBins)+hRawYieldRef->GetBinWidth(nPtBins);
  for(int iPt=0; iPt<nPtBins; iPt++) {
    hMass[iPt] = (TH1F*)infile->Get(Form("hMass_%0.f_%0.f",PtLims[iPt],PtLims[iPt+1]));
  }
  
  //invariant-mass distribution fits
  TH1F* hMassToFit = 0x0;
  TH1F* hRawYield[nPtBins];
  TH1F* hBinCount[nBinCounting][nPtBins];
  TH1F* hSigma[nPtBins];
  TH1F* hMean[nPtBins];
  TH1F* hChiSquare[nPtBins];
  TH1F* hSignal[nPtBins];
  TH1F* hBackground[nPtBins];
  TH1F* hRawYieldVsTrial[nPtBins];
  TH1F* hBinCountVsTrial[nBinCounting][nPtBins];
  TH1F* hSigmaVsTrial[nPtBins];
  TH1F* hMeanVsTrial[nPtBins];
  TH1F* hChiSquareVsTrial[nPtBins];
  
  TNtuple* ntupleMultiTrial = new TNtuple("ntupleMultiTrial","","PtBin:PtMin:PtMax:isSigmaFix:isMeanFix:SgnFunc:BkgFunc:BinWidth:LowFitLimit:HighFitLimit:Raw_yield:Raw_yield_unc:Raw_yield_bincounting1:Raw_yield_bincounting1_unc:Raw_yield_bincounting2:Raw_yield_bincounting2_unc:Sigma:Sigma_unc:Mean:Mean_unc:ChiSquare:Probability:Signal:Signal_unc:Background:Background_unc:Raw_yield_secpeak");
  ntupleMultiTrial->SetDirectory(0);
  
  TLine* lRawRef[nPtBins];
  TLine* lSigmaRef[nPtBins];
  TLine* lMeanRef[nPtBins];
  TLine* lRawRefVsTrial[nPtBins];
  TLine* lSigmaRefVsTrial[nPtBins];
  TLine* lMeanRefVsTrial[nPtBins];
  
  TCanvas* c[nPtBins];
  for(int iPt=0; iPt<nPtBins; iPt++) {
    c[iPt] = new TCanvas(Form("c_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",10,10,1920,1080);
    c[iPt]->Divide(2,2);
  }
  
  double massDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
  double massDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  
  const int nbins=100;
  const int ntrials = nMins*nMaxs*nSigma*nMean*nBkgFunc*nSgnFunc*nReb;

  TPaveText* stats[nPtBins];
  TPaveText* statsbc[nBinCounting][nPtBins];
  TLegend* leg = new TLegend(0.5,0.65,0.89,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  
  int bincountcolors[] = {kOrange+1,kGreen+2,kRed+1};
  int bincountcolorsmarker[] = {kBlack,kCyan,kBlack};
  
  for(int iPt=0; iPt<nPtBins; iPt++) {
    double rawmin = 1.;
    double rawmax = -1.;
    double bkgmin = hBkgRef->GetBinContent(iPt+1+ptshifttoref)*0.5;
    double bkgmax = hBkgRef->GetBinContent(iPt+1+ptshifttoref)*1.5;
    
    if(loadref!=1 && loadref!=2 && loadref!=4) {
      rawmin = hRawYieldRef->GetBinContent(iPt+1+ptshifttoref)*(1-0.6);
      rawmax = hRawYieldRef->GetBinContent(iPt+1+ptshifttoref)*(1+0.6);
      if(PtLims[iPt]>=10) {
        rawmax = hRawYieldRef->GetBinContent(iPt+1+ptshifttoref)*2.;
        rawmin = hRawYieldRef->GetBinContent(iPt+1+ptshifttoref)*0.;
      }
      if(PtLims[iPt]<=3) {
        rawmax = hRawYieldRef->GetBinContent(iPt+1+ptshifttoref)*2.5;
        rawmin = hRawYieldRef->GetBinContent(iPt+1+ptshifttoref)*0.;
      }
      lRawRefVsTrial[iPt] = new TLine(-0.5,hRawYieldRef->GetBinContent(iPt+1+ptshifttoref),ntrials-0.5,hRawYieldRef->GetBinContent(iPt+1+ptshifttoref));
      lRawRefVsTrial[iPt]->SetLineColor(kRed);
      lRawRefVsTrial[iPt]->SetLineWidth(2);
    }
    if(PtLims[iPt]>=36) rawmin=0.;
    double sigmamin = 1.;
    double sigmamax = -1.;
    if(loadref!=1 && loadref!=2 && loadref!=5) {
      sigmamin = hSigmaRef->GetBinContent(iPt+1+ptshifttoref)*(1-0.5);
      sigmamax = hSigmaRef->GetBinContent(iPt+1+ptshifttoref)*(1+0.5);
      if(PtLims[iPt]<3) {
        sigmamin = 0.;
        sigmamax = hSigmaRef->GetBinContent(iPt+1+ptshifttoref)*2.5;
      }
      lSigmaRefVsTrial[iPt] = new TLine(-0.5,hSigmaRef->GetBinContent(iPt+1+ptshifttoref),ntrials-0.5,hSigmaRef->GetBinContent(iPt+1+ptshifttoref));
      lSigmaRefVsTrial[iPt]->SetLineColor(kRed);
      lSigmaRefVsTrial[iPt]->SetLineWidth(2);
    }
    double meanmin = massDs*(1-0.005);
    double meanmax = massDs*(1+0.005);
    if(loadref!=1 && loadref!=2 && loadref!=6) {
      lMeanRefVsTrial[iPt] = new TLine(-0.5,hMeanRef->GetBinContent(iPt+1+ptshifttoref),ntrials-0.5,hMeanRef->GetBinContent(iPt+1+ptshifttoref));
      lMeanRefVsTrial[iPt]->SetLineColor(kRed);
      lMeanRefVsTrial[iPt]->SetLineWidth(2);
    }
    
    hRawYield[iPt] = new TH1F(Form("hRawYield_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",nbins,rawmin,rawmax);
    hSigma[iPt] = new TH1F(Form("hSigma_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",nbins,sigmamin,sigmamax);
    hMean[iPt] = new TH1F(Form("hMean_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",nbins,meanmin,meanmax);
    hChiSquare[iPt] = new TH1F(Form("hChiSquare_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",nbins,0.,maxchisquare);
    hSignal[iPt] = new TH1F(Form("hSignal_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",nbins,rawmin,rawmax);
    hBackground[iPt] = new TH1F(Form("hBackground_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",nbins,bkgmin,bkgmax);
    hRawYield[iPt]->GetXaxis()->SetTitle("raw yield");
    hRawYield[iPt]->GetYaxis()->SetTitle("Entries");
    hRawYield[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hRawYield[iPt]->SetFillStyle(3004);
    hRawYield[iPt]->SetLineWidth(2);
    hRawYield[iPt]->SetLineColor(kBlue+1);
    hRawYield[iPt]->SetFillColor(kBlue+1);
    hSigma[iPt]->GetXaxis()->SetTitle("width (GeV/c^{2})");
    hSigma[iPt]->GetYaxis()->SetTitle("Entries");
    hSigma[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hSigma[iPt]->SetFillStyle(3004);
    hSigma[iPt]->SetLineWidth(2);
    hSigma[iPt]->SetLineColor(kBlue+1);
    hSigma[iPt]->SetFillColor(kBlue+1);
    hMean[iPt]->GetXaxis()->SetTitle("mean (GeV/c^{2})");
    hMean[iPt]->GetYaxis()->SetTitle("Entries");
    hMean[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hMean[iPt]->SetFillStyle(3004);
    hMean[iPt]->SetLineWidth(2);
    hMean[iPt]->SetLineColor(kBlue+1);
    hMean[iPt]->SetFillColor(kBlue+1);
    hChiSquare[iPt]->GetXaxis()->SetTitle("#chi^{2}");
    hChiSquare[iPt]->GetYaxis()->SetTitle("Entries");
    hChiSquare[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hChiSquare[iPt]->SetFillStyle(3004);
    hChiSquare[iPt]->SetLineWidth(2);
    hChiSquare[iPt]->SetLineColor(kBlue+1);
    hChiSquare[iPt]->SetFillColor(kBlue+1);
    hSignal[iPt]->GetXaxis()->SetTitle("signal (3#sigma)");
    hSignal[iPt]->GetYaxis()->SetTitle("Entries");
    hSignal[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hSignal[iPt]->SetFillStyle(3004);
    hSignal[iPt]->SetLineWidth(2);
    hSignal[iPt]->SetLineColor(kBlue+1);
    hSignal[iPt]->SetFillColor(kBlue+1);
    hBackground[iPt]->GetXaxis()->SetTitle("background (3#sigma)");
    hBackground[iPt]->GetYaxis()->SetTitle("Entries");
    hBackground[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hBackground[iPt]->SetFillStyle(3004);
    hBackground[iPt]->SetLineWidth(2);
    hBackground[iPt]->SetLineColor(kBlue+1);
    hBackground[iPt]->SetFillColor(kBlue+1);
    
    hRawYieldVsTrial[iPt] = new TH1F(Form("hRawYieldVsTrial_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",ntrials,-0.5,ntrials-0.5);
    hSigmaVsTrial[iPt] = new TH1F(Form("hMeanVsTrial_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",ntrials,-0.5,ntrials-0.5);
    hMeanVsTrial[iPt] = new TH1F(Form("hSigmaVsTrial_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",ntrials,-0.5,ntrials-0.5);
    hChiSquareVsTrial[iPt] = new TH1F(Form("hChiSquareVsTrial_pT_%0.f-%0.f",PtLims[iPt],PtLims[iPt+1]),"",ntrials,-0.5,ntrials-0.5);
    hRawYieldVsTrial[iPt]->GetYaxis()->SetRangeUser(rawmin-0.2*rawmin,rawmax+0.3*rawmax);
    hRawYieldVsTrial[iPt]->GetYaxis()->SetTitle("raw yield");
    hRawYieldVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hRawYieldVsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hRawYieldVsTrial[iPt]->SetMarkerStyle(20);
    hRawYieldVsTrial[iPt]->SetLineWidth(2);
    hRawYieldVsTrial[iPt]->SetMarkerSize(0.5);
    hRawYieldVsTrial[iPt]->SetLineColor(kBlue+1);
    hRawYieldVsTrial[iPt]->SetMarkerColor(kBlack);
    hSigmaVsTrial[iPt]->GetYaxis()->SetRangeUser(sigmamin,sigmamax);
    hSigmaVsTrial[iPt]->GetYaxis()->SetTitle("width (GeV/c^{2})");
    hSigmaVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hSigmaVsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hSigmaVsTrial[iPt]->SetMarkerStyle(20);
    hSigmaVsTrial[iPt]->SetLineWidth(2);
    hSigmaVsTrial[iPt]->SetMarkerSize(0.5);
    hSigmaVsTrial[iPt]->SetLineColor(kBlue+1);
    hSigmaVsTrial[iPt]->SetMarkerColor(kBlack);
    hMeanVsTrial[iPt]->GetYaxis()->SetRangeUser(meanmin,meanmax);
    hMeanVsTrial[iPt]->GetYaxis()->SetTitle("mean (GeV/c^{2})");
    hMeanVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hMeanVsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hMeanVsTrial[iPt]->SetMarkerStyle(20);
    hMeanVsTrial[iPt]->SetLineWidth(2);
    hMeanVsTrial[iPt]->SetMarkerSize(0.5);
    hMeanVsTrial[iPt]->SetLineColor(kBlue+1);
    hMeanVsTrial[iPt]->SetMarkerColor(kBlack);
    hChiSquareVsTrial[iPt]->GetYaxis()->SetRangeUser(0,maxchisquare);
    hChiSquareVsTrial[iPt]->GetYaxis()->SetTitle("#chi^{2}");
    hChiSquareVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hChiSquareVsTrial[iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
    hChiSquareVsTrial[iPt]->SetMarkerStyle(20);
    hChiSquareVsTrial[iPt]->SetLineWidth(2);
    hChiSquareVsTrial[iPt]->SetMarkerSize(0.5);
    hChiSquareVsTrial[iPt]->SetLineColor(kBlue+1);
    hChiSquareVsTrial[iPt]->SetMarkerColor(kBlack);
    
    for(int iBinCount=0; iBinCount<nBinCounting; iBinCount++) {
      hBinCount[iBinCount][iPt] = new TH1F(Form("hBinCount_%0.1fsigma_pT_%0.f-%0.f",nSigmaBinCount[iBinCount],PtLims[iPt],PtLims[iPt+1]),"",nbins,rawmin,rawmax);
      hBinCount[iBinCount][iPt]->GetXaxis()->SetTitle("raw yield");
      hBinCount[iBinCount][iPt]->GetYaxis()->SetTitle("Entries");
      hBinCount[iBinCount][iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
      hBinCount[iBinCount][iPt]->SetFillStyle(3004);
      hBinCount[iBinCount][iPt]->SetLineWidth(2);
      hBinCount[iBinCount][iPt]->SetLineColor(bincountcolors[iBinCount]);
      hBinCount[iBinCount][iPt]->SetFillColor(bincountcolors[iBinCount]);
      
      hBinCountVsTrial[iBinCount][iPt] = new TH1F(Form("hBinCountVsTrial_%0.1fsigma_pT_%0.f-%0.f",nSigmaBinCount[iBinCount],PtLims[iPt],PtLims[iPt+1]),"",ntrials,-0.5,ntrials-0.5);
      {hBinCountVsTrial[iBinCount][iPt]->GetYaxis()->SetRangeUser(rawmin-0.2*rawmin,rawmax+0.3*rawmax);}
      hBinCountVsTrial[iBinCount][iPt]->GetYaxis()->SetTitle("raw yield");
      hBinCountVsTrial[iBinCount][iPt]->GetXaxis()->SetTitle("Trial #");
      hBinCountVsTrial[iBinCount][iPt]->SetTitle(Form("%0.f < #it{p}_{T} < %0.f (GeV/c)",PtLims[iPt],PtLims[iPt+1]));
      hBinCountVsTrial[iBinCount][iPt]->SetMarkerStyle(20);
      hBinCountVsTrial[iBinCount][iPt]->SetLineWidth(2);
      hBinCountVsTrial[iBinCount][iPt]->SetMarkerSize(0.5);
      hBinCountVsTrial[iBinCount][iPt]->SetLineColor(bincountcolors[iBinCount]);
      hBinCountVsTrial[iBinCount][iPt]->SetMarkerColor(bincountcolorsmarker[iBinCount]);
    }
        
    int iTrial=0;
    
    for(int iSgnFunc=0; iSgnFunc<nSgnFunc; iSgnFunc++) {
      for(int iSigma=kFree; iSigma<nSigma; iSigma++) {
        for(int iMean=kFree; iMean<nMean; iMean++) {
          for(int iBkgFunc=0; iBkgFunc<nBkgFunc; iBkgFunc++) {
            for(int iReb=0; iReb<nReb; iReb++) {
              for(int iMin=0; iMin<nMins; iMin++) {
                for(int iMax=0; iMax<nMaxs; iMax++) {
                  
                  hMassToFit=(TH1F*)AliVertexingHFUtils::RebinHisto(hMass[iPt],rebin[iReb]);     
                  AliHFInvMassFitter fitter(hMassToFit,mins[iMin],maxs[iMax],bkgfcn[iBkgFunc],sgnfcn[iSgnFunc]);
                  fitter.SetUseLikelihoodFit();
                  
                  if(iMean==kFixed) fitter.SetFixGaussianMean(massDs);
                  else fitter.SetInitialGaussianMean(massDs);
                  
                  if(sgnfcn[iSgnFunc]==AliHFInvMassFitter::k2GausSigmaRatioPar) {
                    if(loadref>8) continue;
                    double frac2gaus = hFracGaus2MC->GetBinContent(iPt+1);
                    double ratio2gaussigma = hRatioSigmaGaus2MC->GetBinContent(iPt+1);
                    fitter.SetFixFrac2Gaus(frac2gaus);
                    fitter.SetFixRatio2GausSigma(ratio2gaussigma);
                  }

                  double MCsigma=-1;
                  if(loadref!=2 && loadref!=3 && loadref!=7) {
                    MCsigma = hSigmaMC->GetBinContent(iPt+1+ptshifttoref);
                    if(iSigma==kFixed) fitter.SetFixGaussianSigma(MCsigma);
                    else if(iSigma==kFixedPlus15Perc) {MCsigma += 0.15*MCsigma; fitter.SetFixGaussianSigma(MCsigma);}
                    else if(iSigma==kFixedMinus15Perc) {MCsigma -= 0.15*MCsigma; fitter.SetFixGaussianSigma(MCsigma);}
                    else fitter.SetInitialGaussianSigma(MCsigma);
                  }
                  else {
                    if(iSigma==kFixed) continue;
                    else fitter.SetInitialGaussianSigma(0.011);
                  }
                  fitter.IncludeSecondGausPeak(massDplus,kTRUE,0.008,kTRUE);
                  int fitok = fitter.MassFitter(kFALSE);
                  if(fitok==0) continue; //failed fit
                  
                  double chi=fitter.GetReducedChiSquare();
                  double prob=fitter.GetFitProbability();
                  double rawyield=fitter.GetRawYield();
                  double sigma=fitter.GetSigma();
                  double mean=fitter.GetMean();
                  double rawyielderr=fitter.GetRawYieldError();
                  double sigmaerr=fitter.GetSigmaUncertainty();
                  double meanerr=fitter.GetMeanUncertainty();
                  
                  double bkg=-1, bkgerr=-1;
                  double sgn=-1, sgnerr=-1;
                  fitter.Background(3,bkg,bkgerr);
                  fitter.Signal(3,sgn,sgnerr);
                  
                  double cntSig[nBinCounting];
                  double cntErr[nBinCounting];
                  for(int iBinCount=0; iBinCount<nBinCounting; iBinCount++) cntSig[iBinCount]=fitter.GetRawYieldBinCounting(cntErr[iBinCount], nSigmaBinCount[iBinCount],1);
                
                  TF1* fTotFunc = (TF1*)fitter.GetMassFunc();
                  double rawyieldsecpeak = fTotFunc->GetParameter(fTotFunc->GetNpar()-3)/hMassToFit->GetBinWidth(1);

                  //only 2 bin counting stored in ntuple
                  float arrayforntuple[27] = {(float)iPt+ptshifttoref,(float)PtLims[iPt],(float)PtLims[iPt+1],(float)iSigma,(float)iMean,(float)sgnfcn[iSgnFunc],(float)bkgfcn[iBkgFunc],(float)hMassToFit->GetBinWidth(1),(float)mins[iMin],(float)maxs[iMax],(float)rawyield,(float)rawyielderr,(float)cntSig[0],(float)cntErr[0],(float)cntSig[1],(float)cntErr[1],(float)sigma,(float)sigmaerr,(float)mean,(float)meanerr,(float)chi,(float)prob,(float)sgn,(float)sgnerr,(float)bkg,(float)bkgerr,(float)rawyieldsecpeak};
                  ntupleMultiTrial->Fill(arrayforntuple);

                  if(chi<maxchisquare && chi>minchisquare) {
                    hRawYield[iPt]->Fill(rawyield);
                    hSigma[iPt]->Fill(sigma);
                    hMean[iPt]->Fill(mean);
                    hChiSquare[iPt]->Fill(chi);
                    hSignal[iPt]->Fill(sgn);
                    hBackground[iPt]->Fill(bkg);
                    hRawYieldVsTrial[iPt]->SetBinContent(iTrial+1,rawyield);
                    hSigmaVsTrial[iPt]->SetBinContent(iTrial+1,sigma);
                    hMeanVsTrial[iPt]->SetBinContent(iTrial+1,mean);
                    hChiSquareVsTrial[iPt]->SetBinContent(iTrial+1,chi);
                    hRawYieldVsTrial[iPt]->SetBinError(iTrial+1,rawyielderr);
                    hSigmaVsTrial[iPt]->SetBinError(iTrial+1,sigmaerr);
                    hMeanVsTrial[iPt]->SetBinError(iTrial+1,meanerr);
                    
                    for(int iBinCount=0; iBinCount<nBinCounting; iBinCount++) {
                      hBinCount[iBinCount][iPt]->Fill(cntSig[iBinCount]);
                      hBinCountVsTrial[iBinCount][iPt]->SetBinContent(iTrial+1,cntSig[iBinCount]);
                      hBinCountVsTrial[iBinCount][iPt]->SetBinError(iTrial+1,cntErr[iBinCount]);
                    }
                  }
                  
                  iTrial++;
                  
                  delete hMassToFit;
                  hMassToFit=0x0;
                }
              }
            }
          }
        }
      }
    }
    
    double mean = hRawYield[iPt]->GetMean();
    double RMS = hRawYield[iPt]->GetRMS();
    double disp = (max(hRawYieldVsTrial[iPt])-min(hRawYieldVsTrial[iPt]))/TMath::Sqrt(12);
    double RMSunc = RMS/hRawYieldRef->GetBinContent(iPt+1+ptshifttoref)*100;
    double Flatunc = disp/hRawYieldRef->GetBinContent(iPt+1+ptshifttoref)*100;
    
    stats[iPt]=new TPaveText(0.15,0.65,0.44,0.85,"NDC");
    stats[iPt]->SetTextSize(0.05);
    stats[iPt]->SetFillColor(0);
    stats[iPt]->SetFillStyle(0);
    stats[iPt]->SetBorderSize(0);
    stats[iPt]->SetTextFont(42);
    stats[iPt]->SetTextColor(kBlue+1);
    stats[iPt]->AddText(Form("mean = %0.1f",mean));
    stats[iPt]->AddText(Form("RMS = %0.1f (%0.1f%%)",RMS,RMSunc));
    stats[iPt]->AddText(Form("#frac{max-min}{#sqrt{12}} = %.1f (%0.1f%%)",disp,Flatunc));
    
    for(int iBinCount=0; iBinCount<nBinCounting; iBinCount++) {
      
      mean = hBinCount[iBinCount][iPt]->GetMean();
      RMS = hBinCount[iBinCount][iPt]->GetRMS();
      disp = (max(hBinCountVsTrial[iBinCount][iPt])-min(hBinCountVsTrial[iBinCount][iPt]))/TMath::Sqrt(12);
      RMSunc = RMS/hRawYieldRef->GetBinContent(iPt+1+ptshifttoref)*100;
      Flatunc = disp/hRawYieldRef->GetBinContent(iPt+1+ptshifttoref)*100;
      
      statsbc[iBinCount][iPt]=new TPaveText(0.6,0.65-0.3*iBinCount,0.89,0.85-0.3*iBinCount,"NDC");
      statsbc[iBinCount][iPt]->SetTextSize(0.05);
      statsbc[iBinCount][iPt]->SetFillColor(0);
      statsbc[iBinCount][iPt]->SetFillStyle(0);
      statsbc[iBinCount][iPt]->SetBorderSize(0);
      statsbc[iBinCount][iPt]->SetTextFont(42);
      statsbc[iBinCount][iPt]->SetTextColor(bincountcolors[iBinCount]);
      statsbc[iBinCount][iPt]->AddText(Form("mean = %0.1f",mean));
      statsbc[iBinCount][iPt]->AddText(Form("RMS = %0.1f (%0.1f%%)",RMS, RMSunc));
      statsbc[iBinCount][iPt]->AddText(Form("#frac{max-min}{#sqrt{12}} = %0.1f (%0.1f%%)",disp,Flatunc));
    }
    
    if(iPt==0) {
      leg->AddEntry(lRawRefVsTrial[iPt],"Central value","l");
      leg->AddEntry(hRawYieldVsTrial[iPt],"Fit method","lpe");
      for(int iBinCount=0; iBinCount<nBinCounting; iBinCount++) {
        leg->AddEntry(hBinCountVsTrial[iBinCount][iPt],Form("Bin counting (%0.1f#sigma)",nSigmaBinCount[iBinCount]),"lpe");
      }
    }
    
    c[iPt]->cd(1)->SetTopMargin(0.12);
    hBinCountVsTrial[0][iPt]->Draw();
    for(int iBinCount=1; iBinCount<nBinCounting; iBinCount++) {
      hBinCountVsTrial[iBinCount][iPt]->Draw("same");
    }
    hRawYieldVsTrial[iPt]->Draw("same");
    if(lRawRefVsTrial[iPt]) {lRawRefVsTrial[iPt]->Draw("same");}
    leg->Draw("same");
    c[iPt]->cd(2)->SetTopMargin(0.12);;
    hBinCount[0][iPt]->GetYaxis()->SetRangeUser(0.,hRawYield[iPt]->GetMaximum()*1.5);
    hBinCount[0][iPt]->Draw();
    for(int iBinCount=1; iBinCount<nBinCounting; iBinCount++) {
      hBinCount[iBinCount][iPt]->Draw("same");
    }
    hRawYield[iPt]->Draw("same");
    if(loadref!=1 && loadref!=2 && loadref!=4) {
      lRawRef[iPt] = new TLine(hRawYieldRef->GetBinContent(iPt+1+ptshifttoref),0,hRawYieldRef->GetBinContent(iPt+1+ptshifttoref),hRawYield[iPt]->GetMaximum()*1.5);
      lRawRef[iPt]->SetLineColor(kRed);
      lRawRef[iPt]->SetLineWidth(2);
      lRawRef[iPt]->Draw("same");
    }
    stats[iPt]->Draw("same");
    for(int iBinCount=0; iBinCount<nBinCounting; iBinCount++) {
      statsbc[iBinCount][iPt]->Draw("same");
    }
    c[iPt]->cd(3)->SetTopMargin(0.12);;
    hSigmaVsTrial[iPt]->Draw();
    if(lSigmaRefVsTrial[iPt]) {lSigmaRefVsTrial[iPt]->Draw("same");}
    c[iPt]->cd(4)->SetTopMargin(0.12);;
    hChiSquareVsTrial[iPt]->Draw("p");
  }
  
  //output files
  TFile outfileraw(outfilerawname.Data(),"RECREATE");
  ntupleMultiTrial->Write();
  for(int iPt=0; iPt<nPtBins; iPt++) {
    hRawYield[iPt]->Write();
    hSigma[iPt]->Write();
    hMean[iPt]->Write();
    hChiSquare[iPt]->Write();
    hSignal[iPt]->Write();
    hBackground[iPt]->Write();
    hRawYieldVsTrial[iPt]->Write();
    hSigmaVsTrial[iPt]->Write();
    hMeanVsTrial[iPt]->Write();
    hChiSquareVsTrial[iPt]->Write();
    for(int iBinCount=0; iBinCount<nBinCounting; iBinCount++) {
      hBinCount[iBinCount][iPt]->Write();
      hBinCountVsTrial[iBinCount][iPt]->Write();
    }
    c[iPt]->Write();
  }
  outfileraw.Close();
  cout << "\n" <<outfilerawname << " saved." <<endl;
  
  for(int iPt=0; iPt<nPtBins; iPt++) {
    outfilerawname.ReplaceAll(".root","");
    c[iPt]->SaveAs(Form("%s_pt%d.pdf",outfilerawname.Data(),iPt+ptshifttoref));
  }
    
  return 0;
}

//__________________________________________________________________________________________________________________
double min(TH1F* histo) {
  
  double min=1.e+20;
  for(int iBin=0; iBin<histo->GetNbinsX(); iBin++) {
    if(histo->GetBinContent(iBin+1)>0 && histo->GetBinContent(iBin+1)<min) {min = histo->GetBinContent(iBin+1);}
  }
  return min;
}

//__________________________________________________________________________________________________________________
double max(TH1F* histo) {
  
  double max=-1.e+20;
  for(int iBin=0; iBin<histo->GetNbinsX(); iBin++) {
    if(histo->GetBinContent(iBin+1)>0 && histo->GetBinContent(iBin+1)>max) {max = histo->GetBinContent(iBin+1);}
  }
  return max;
}

//__________________________________________________________________________________________________________________
int LoadRefFiles(TString reffilename, TString reffileMCname, TString reffileMCname2, TH1F *&hRawYieldRef, TH1F *&hSigmaRef, TH1F *&hMeanRef, TH1F *&hBkgRef, TH1F *&hSigmaMC, TH1F *&hFracGaus2MC, TH1F *&hSigmaGaus2MC) {

  TFile* reffile = TFile::Open(reffilename.Data(),"READ");
  if(reffile) {
    hRawYieldRef=(TH1F*)reffile->Get("hRawYields");
    hSigmaRef=(TH1F*)reffile->Get("hRawYieldsSigma");
    hMeanRef=(TH1F*)reffile->Get("hRawYieldsMean");
    hBkgRef=(TH1F*)reffile->Get("hRawYieldsBkg");
    if(hRawYieldRef) {hRawYieldRef->SetDirectory(0);}
    if(hSigmaRef) {hSigmaRef->SetDirectory(0);}
    if(hMeanRef) {hMeanRef->SetDirectory(0);}
    if(hBkgRef) {hBkgRef->SetDirectory(0);}
    reffile->Close();
  }

  TFile* reffileMC = TFile::Open(reffileMCname.Data(),"READ");
  if(reffileMC) {
    hSigmaMC=(TH1F*)reffileMC->Get("hRawYieldsSigma");
    if(hSigmaMC) {hSigmaMC->SetDirectory(0);}
    reffileMC->Close();
  }

  TFile* reffileMC2 = TFile::Open(reffileMCname2.Data(),"READ");
  if(reffileMC2) {
    hFracGaus2MC=(TH1F*)reffileMC2->Get("hRawYieldsFracGaus2");
    if(hFracGaus2MC) {hFracGaus2MC->SetDirectory(0);}
    hSigmaGaus2MC=(TH1F*)reffileMC2->Get("hRawYieldsSigma2");
    if(hSigmaGaus2MC) {hSigmaGaus2MC->SetDirectory(0);}
    reffileMC2->Close();
  }

  if(!reffile && !reffileMC) return 1;
  if(!reffile) return 2;
  if(!reffileMC) return 3;
  if(!hRawYieldRef) return 4;
  if(!hSigmaRef) return 5;
  if(!hMeanRef) return 6;
  if(!hBkgRef) return 7;
  if(!hSigmaMC) return 8;
  if(!hFracGaus2MC) return 9;
  if(!hSigmaGaus2MC) return 10;
  return 0;
}

//__________________________________________________________________________________________________________________
void SetStyle() {
  gStyle->SetOptFit(1);
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
