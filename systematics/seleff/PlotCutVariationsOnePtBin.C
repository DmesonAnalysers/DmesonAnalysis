#if !defined (__CINT__) || defined (__CLING__)

#include <Riostream.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <THnSparse.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TAttMarker.h>
#include <TGaxis.h>
#include <TGraphErrors.h>

#endif

const TString inputdirname = "../ML_DsAnalysis/QM_prel/outputs/3050/eff_syst_highpt";
const TString inputfilecommonname_crosssec = "DsYield_method2_fd1_br1";
const TString inputfilecommonname_rawyield =  "RawYieldsDs";
const TString inputfilecommonname_eff = "Efficiency_Ds";
const TString inputfilesuffix[] = {
                                   };

double maxchi2 = 2.;
double minsignif = 3;
double minrelsignif = 0.5;
double minreleff = 0.3;
double maxreleff = 1.65;
double fillthrreleff = 0.005;

const  double relassignedsyst[] = {0.06,0.06,0.06,0.05};

const TString outfilename = "../ML_DsAnalysis/QM_prel/outputs/3050/eff_syst_lowpt/CutVarSyst_Ds_3050_hight.root";

int PlotCutVariationsOnePtBin( bool fRelativeVariation=kFALSE);
void SetStyle();

 int PlotCutVariationsOnePtBin( bool fRelativeVariation) {
  
  SetStyle();
  gStyle->SetTitleSize(0.05,"xy");
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadRightMargin(0.12);

  const  int nFiles = sizeof(inputfilesuffix)/sizeof(inputfilesuffix[0]);
  TH1D* hCrossSection[nFiles];
  TH1D* hRawYield[nFiles];
  TH1D* hSignificance[nFiles];
  TH1D* hSoverB[nFiles];
  TH1D* hchi2[nFiles];
  TH1D* hEffPrompt[nFiles];
  TH1D* hEffFD[nFiles];

  for( int iFile=0; iFile<nFiles; iFile++) {
    TFile* infile_crossec = TFile::Open(Form("%s/%s%s.root",inputdirname.Data(),Form("crosssec/%s",inputfilecommonname_crosssec.Data()),inputfilesuffix[iFile].Data()));
    if(!infile_crossec) return iFile+1;
    hCrossSection[iFile] = (TH1D*)infile_crossec->Get("hAAC");
    hCrossSection[iFile]->SetDirectory(0);
    TFile* infile_rawyield = TFile::Open(Form("%s/%s%s.root",inputdirname.Data(),Form("rawyields/%s",inputfilecommonname_rawyield.Data()),inputfilesuffix[iFile].Data()));
    if(!infile_rawyield) return iFile+11;
    hRawYield[iFile] = (TH1D*)infile_rawyield->Get("hRawYields");
    hSignificance[iFile] = (TH1D*)infile_rawyield->Get("hRawYieldsSignificance");
    hSoverB[iFile] = (TH1D*)infile_rawyield->Get("hRawYieldsSoverB");
    hchi2[iFile] = (TH1D*)infile_rawyield->Get("hRawYieldsChiSquare");
    hRawYield[iFile]->SetDirectory(0);
    hSignificance[iFile]->SetDirectory(0);
    hSoverB[iFile]->SetDirectory(0);
    hchi2[iFile]->SetDirectory(0);
    infile_rawyield->Close();
    TFile* infile_eff = TFile::Open(Form("%s/%s%s.root",inputdirname.Data(),Form("efficiency/%s",inputfilecommonname_eff.Data()),inputfilesuffix[iFile].Data()));
    if(!infile_eff) return iFile+111;
    hEffPrompt[iFile] = (TH1D*)infile_eff->Get("hEffPrompt");
    hEffFD[iFile] = (TH1D*)infile_eff->Get("hEffFD");
    hEffPrompt[iFile]->SetDirectory(0);
    hEffFD[iFile]->SetDirectory(0);
    infile_eff->Close();
  }
  
  const  int nPtBins = hRawYield[0]->GetNbinsX();
  TGraphErrors* gCrossSectionVsCutSet[nPtBins];
  TGraphErrors* gRawYieldVsCutSet[nPtBins];
  TGraphErrors* gSignificanceVsCutSet[nPtBins];
  TGraphErrors* gSoverBVsCutSet[nPtBins];
  TGraphErrors* gEffPromptVsCutSet[nPtBins];
  TGraphErrors* gEffFDVsCutSet[nPtBins];
  TGraphErrors* gCrossSectionCent[nPtBins];
  TH1D* hCrossSectionRatioDist[nPtBins];
  TLine* lCrossSectionCent[nPtBins];
  TBox* bCrossSectionSyst[nPtBins];
  TBox* bCrossSectionRMS[nPtBins];
  TCanvas* cOutPut[nPtBins];

   double maxrawyield[nPtBins];
   double maxefficiency[nPtBins];
   double maxsignif[nPtBins];
   double maxSoverB[nPtBins];
   double maxCross[nPtBins];

  TLegend* legEff = new TLegend(0.2,0.7,0.5,0.85);
  legEff->SetTextSize(0.05);
  legEff->SetFillStyle(0);
  legEff->SetBorderSize(0);

  TLegend* legCross = new TLegend(0.2,0.75,0.45,0.85);
  legCross->SetTextSize(0.05);
  legCross->SetFillStyle(0);
  legCross->SetBorderSize(0);

  TLegend* legSyst = new TLegend(0.2,0.7,0.45,0.85);
  legSyst->SetTextSize(0.05);
  legSyst->SetFillStyle(0);
  legSyst->SetBorderSize(0);

  for( int iPtRaw=0; iPtRaw<nPtBins; iPtRaw++) {
    cOutPut[iPtRaw] = new TCanvas(Form("cOutPut_ptbin%d",iPtRaw),"",1920,1080);
    cOutPut[iPtRaw]->Divide(3,2);

    gCrossSectionCent[iPtRaw] = new TGraphErrors(1);
    gCrossSectionCent[iPtRaw]->SetTitle(";Cut set; d#it{N}/d#it{p}_{T} (GeV^{-1} #it{c})");
    gCrossSectionCent[iPtRaw]->SetName(Form("gCrossSectionCent_ptbin%d",iPtRaw));
    gCrossSectionCent[iPtRaw]->SetLineColor(kRed+1);
    gCrossSectionCent[iPtRaw]->SetFillColorAlpha(kRed+1,0.2);
    gCrossSectionCent[iPtRaw]->SetLineWidth(2);
    
    gCrossSectionVsCutSet[iPtRaw] = new TGraphErrors(0);
    gCrossSectionVsCutSet[iPtRaw]->SetTitle(";Cut set; d#it{N}/d#it{p}_{T} (GeV^{-1} #it{c})");
    gCrossSectionVsCutSet[iPtRaw]->SetName(Form("gCrossSectionVsCutSet_ptbin%d",iPtRaw));
    gCrossSectionVsCutSet[iPtRaw]->SetMarkerSize(1.);
    gCrossSectionVsCutSet[iPtRaw]->SetMarkerStyle(kFullCircle);
    gCrossSectionVsCutSet[iPtRaw]->SetMarkerColor(kBlack);
    gCrossSectionVsCutSet[iPtRaw]->SetLineColor(kBlack);
    gCrossSectionVsCutSet[iPtRaw]->SetLineWidth(2);

    gRawYieldVsCutSet[iPtRaw] = new TGraphErrors(0);
    gRawYieldVsCutSet[iPtRaw]->SetTitle(";Cut set; raw yield");
    gRawYieldVsCutSet[iPtRaw]->SetName(Form("gRawYieldVsCutSet_ptbin%d",iPtRaw));
    gRawYieldVsCutSet[iPtRaw]->SetMarkerSize(1.);
    gRawYieldVsCutSet[iPtRaw]->SetMarkerStyle(kFullCircle);
    gRawYieldVsCutSet[iPtRaw]->SetMarkerColor(kBlack);
    gRawYieldVsCutSet[iPtRaw]->SetLineColor(kBlack);
    gRawYieldVsCutSet[iPtRaw]->SetLineWidth(2);

    gSignificanceVsCutSet[iPtRaw] = new TGraphErrors(0);
    gSignificanceVsCutSet[iPtRaw]->SetTitle(";Cut set; significance");
    gSignificanceVsCutSet[iPtRaw]->SetName(Form("gSignificanceVsCutSet_ptbin%d",iPtRaw));
    gSignificanceVsCutSet[iPtRaw]->SetMarkerSize(1.);
    gSignificanceVsCutSet[iPtRaw]->SetMarkerStyle(kFullCircle);
    gSignificanceVsCutSet[iPtRaw]->SetMarkerColor(kBlack);
    gSignificanceVsCutSet[iPtRaw]->SetLineColor(kBlack);
    gSignificanceVsCutSet[iPtRaw]->SetLineWidth(2);

    gSoverBVsCutSet[iPtRaw] = new TGraphErrors(0);
    gSoverBVsCutSet[iPtRaw]->SetTitle(";Cut set; S/B (3#sigma)");
    gSoverBVsCutSet[iPtRaw]->SetName(Form("gSoverBVsCutSet_ptbin%d",iPtRaw));
    gSoverBVsCutSet[iPtRaw]->SetMarkerSize(1.);
    gSoverBVsCutSet[iPtRaw]->SetMarkerStyle(kFullCircle);
    gSoverBVsCutSet[iPtRaw]->SetMarkerColor(kBlack);
    gSoverBVsCutSet[iPtRaw]->SetLineColor(kBlack);
    gSoverBVsCutSet[iPtRaw]->SetLineWidth(2);

    gEffPromptVsCutSet[iPtRaw] = new TGraphErrors(0);
    gEffPromptVsCutSet[iPtRaw]->SetTitle(";Cut set; Efficiency");
    gEffPromptVsCutSet[iPtRaw]->SetName(Form("gEffPromptVsCutSet_ptbin%d",iPtRaw));
    gEffPromptVsCutSet[iPtRaw]->SetMarkerSize(1.);
    gEffPromptVsCutSet[iPtRaw]->SetMarkerStyle(kFullCircle);
    gEffPromptVsCutSet[iPtRaw]->SetMarkerColor(kRed);
    gEffPromptVsCutSet[iPtRaw]->SetLineColor(kRed);
    gEffPromptVsCutSet[iPtRaw]->SetLineWidth(2);

    gEffFDVsCutSet[iPtRaw] = new TGraphErrors(0);
    gEffFDVsCutSet[iPtRaw]->SetTitle(";Cut set; Efficiency");
    gEffFDVsCutSet[iPtRaw]->SetName(Form("gEffFDVsCutSet_ptbin%d",iPtRaw));
    gEffFDVsCutSet[iPtRaw]->SetMarkerSize(1.);
    gEffFDVsCutSet[iPtRaw]->SetMarkerStyle(kFullSquare);
    gEffFDVsCutSet[iPtRaw]->SetMarkerColor(kBlue);
    gEffFDVsCutSet[iPtRaw]->SetLineColor(kBlue);
    gEffFDVsCutSet[iPtRaw]->SetLineWidth(2);

    hCrossSectionRatioDist[iPtRaw] = new TH1D(Form("hCrossSectionRatioDist_ptbin%d",iPtRaw),";(d#it{N}/d#it{p}_{T}) / (d#it{N}/d#it{p}_{T})_{central};Entries",60,0.45,1.55);
    hCrossSectionRatioDist[iPtRaw]->SetLineColor(kBlack);
    hCrossSectionRatioDist[iPtRaw]->SetLineWidth(2);

    maxrawyield[iPtRaw]=-1;
    maxefficiency[iPtRaw]=-1;
    maxsignif[iPtRaw]=-1;
    maxSoverB[iPtRaw]=-1;
    maxCross[iPtRaw]=-1;

    for( int iFile=0; iFile<nFiles; iFile++) {
      if(maxrawyield[iPtRaw]<hRawYield[iFile]->GetBinContent(iPtRaw+1)) maxrawyield[iPtRaw] = hRawYield[iFile]->GetBinContent(iPtRaw+1);
      if(maxefficiency[iPtRaw]<hEffFD[iFile]->GetBinContent(iPtRaw+1)) maxefficiency[iPtRaw] = hEffFD[iFile]->GetBinContent(iPtRaw+1);
      if(maxsignif[iPtRaw]<hSignificance[iFile]->GetBinContent(iPtRaw+1)) maxsignif[iPtRaw] = hSignificance[iFile]->GetBinContent(iPtRaw+1);
      if(maxSoverB[iPtRaw]<hSoverB[iFile]->GetBinContent(iPtRaw+1)) maxSoverB[iPtRaw] = hSoverB[iFile]->GetBinContent(iPtRaw+1);
      if(maxCross[iPtRaw]<hCrossSection[iFile]->GetBinContent(iPtRaw+1)) maxCross[iPtRaw] = hCrossSection[iFile]->GetBinContent(iPtRaw+1);
      
       int crossptbin = hCrossSection[iFile]->GetXaxis()->FindBin(hRawYield[iFile]->GetBinCenter(iPtRaw+1));
       int effptbin = hEffPrompt[iFile]->GetXaxis()->FindBin(hRawYield[iFile]->GetBinCenter(iPtRaw+1));

      if(iFile==0) {
        lCrossSectionCent[iPtRaw] = new TLine(0.,hCrossSection[iFile]->GetBinContent(crossptbin),nFiles-1,hCrossSection[iFile]->GetBinContent(crossptbin));
        lCrossSectionCent[iPtRaw]->SetLineWidth(2);
        lCrossSectionCent[iPtRaw]->SetLineColor(kRed+1);

        gCrossSectionCent[iPtRaw]->SetPoint(iFile,nFiles/2,hCrossSection[iFile]->GetBinContent(crossptbin));
        gCrossSectionCent[iPtRaw]->SetPointError(iFile,nFiles/2,hCrossSection[iFile]->GetBinContent(crossptbin)*relassignedsyst[iPtRaw]);
      }
          
      gCrossSectionVsCutSet[iPtRaw]->SetPoint(iFile,iFile,hCrossSection[iFile]->GetBinContent(crossptbin));
      gCrossSectionVsCutSet[iPtRaw]->SetPointError(iFile,0.5,hCrossSection[iFile]->GetBinError(crossptbin));

      if(!fRelativeVariation) {
       gRawYieldVsCutSet[iPtRaw]->SetPoint(iFile,iFile,hRawYield[iFile]->GetBinContent(iPtRaw+1));
       gRawYieldVsCutSet[iPtRaw]->SetPointError(iFile,0.5,hRawYield[iFile]->GetBinError(iPtRaw+1));

       gEffPromptVsCutSet[iPtRaw]->SetPoint(iFile,iFile,hEffPrompt[iFile]->GetBinContent(effptbin));
       gEffPromptVsCutSet[iPtRaw]->SetPointError(iFile,0.5,hEffPrompt[iFile]->GetBinError(effptbin));

       gEffFDVsCutSet[iPtRaw]->SetPoint(iFile,iFile,hEffFD[iFile]->GetBinContent(effptbin));
       gEffFDVsCutSet[iPtRaw]->SetPointError(iFile,0.5,hEffFD[iFile]->GetBinError(effptbin));
      }
      else {
       gRawYieldVsCutSet[iPtRaw]->SetPoint(iFile,iFile,hRawYield[iFile]->GetBinContent(iPtRaw+1)/hRawYield[0]->GetBinContent(iPtRaw+1));
       gRawYieldVsCutSet[iPtRaw]->SetPointError(iFile,0.5,0);

       gEffPromptVsCutSet[iPtRaw]->SetPoint(iFile,iFile,hEffPrompt[iFile]->GetBinContent(effptbin)/hEffPrompt[0]->GetBinContent(effptbin));
       gEffPromptVsCutSet[iPtRaw]->SetPointError(iFile,0.5,0);

       gEffFDVsCutSet[iPtRaw]->SetPoint(iFile,iFile,hEffFD[iFile]->GetBinContent(effptbin)/hEffFD[0]->GetBinContent(effptbin));
       gEffFDVsCutSet[iPtRaw]->SetPointError(iFile,0.5,0);
      }

      gSignificanceVsCutSet[iPtRaw]->SetPoint(iFile,iFile,hSignificance[iFile]->GetBinContent(iPtRaw+1));
      gSignificanceVsCutSet[iPtRaw]->SetPointError(iFile,0.5,hSignificance[iFile]->GetBinError(iPtRaw+1));

      gSoverBVsCutSet[iPtRaw]->SetPoint(iFile,iFile,hSoverB[iFile]->GetBinContent(iPtRaw+1));
      gSoverBVsCutSet[iPtRaw]->SetPointError(iFile,0.5,hSoverB[iFile]->GetBinError(iPtRaw+1));      
    
      if(iFile!=0) {
        if(hchi2[iFile]->GetBinContent(iPtRaw+1) > maxchi2) continue;
        if(TMath::Abs(1-hEffPrompt[iFile]->GetBinContent(effptbin)/hEffPrompt[0]->GetBinContent(effptbin)) < fillthrreleff && 
           TMath::Abs(1-hEffFD[iFile]->GetBinContent(iPtRaw+1)/hEffFD[0]->GetBinContent(iPtRaw+1)) < fillthrreleff) continue;
        if((hEffPrompt[iFile]->GetBinContent(effptbin)/hEffPrompt[0]->GetBinContent(effptbin)) < minreleff) continue;
        if((hEffPrompt[iFile]->GetBinContent(effptbin)/hEffPrompt[0]->GetBinContent(effptbin)) > maxreleff) continue;
        if(hSignificance[iFile]->GetBinContent(iPtRaw+1) < minsignif) continue;
        if((hSignificance[iFile]->GetBinContent(iPtRaw+1)/hSignificance[0]->GetBinContent(iPtRaw+1)) < minrelsignif) continue;
        hCrossSectionRatioDist[iPtRaw]->Fill(hCrossSection[iFile]->GetBinContent(crossptbin)/hCrossSection[0]->GetBinContent(crossptbin));
      }
    }

    bCrossSectionSyst[iPtRaw] = new TBox(1-relassignedsyst[iPtRaw],0.,1+relassignedsyst[iPtRaw],hCrossSectionRatioDist[iPtRaw]->GetMaximum());
    bCrossSectionSyst[iPtRaw]->SetLineColor(kRed+1);
    bCrossSectionSyst[iPtRaw]->SetFillColorAlpha(kRed+1,0.2);

    bCrossSectionRMS[iPtRaw] = new TBox(1-TMath::Sqrt(hCrossSectionRatioDist[iPtRaw]->GetRMS()*hCrossSectionRatioDist[iPtRaw]->GetRMS()+(1-hCrossSectionRatioDist[iPtRaw]->GetMean())*(1-hCrossSectionRatioDist[iPtRaw]->GetMean())),0.,1+TMath::Sqrt(hCrossSectionRatioDist[iPtRaw]->GetRMS()*hCrossSectionRatioDist[iPtRaw]->GetRMS()+(1-hCrossSectionRatioDist[iPtRaw]->GetMean())*(1-hCrossSectionRatioDist[iPtRaw]->GetMean())),hCrossSectionRatioDist[iPtRaw]->GetMaximum());
    bCrossSectionRMS[iPtRaw]->SetLineColor(kBlue+1);
    bCrossSectionRMS[iPtRaw]->SetFillColorAlpha(kBlue+1,0.2);

    if(iPtRaw==0) {
      legEff->AddEntry(gEffPromptVsCutSet[iPtRaw],"Prompt","lpe");
      legEff->AddEntry(gEffFDVsCutSet[iPtRaw],"Feed-down","lpe");
      legCross->AddEntry(gCrossSectionCent[iPtRaw],"Central value #pm assigned syst","fl");
      legSyst->AddEntry(bCrossSectionSyst[iPtRaw],"Assigned syst","f");
      legSyst->AddEntry(bCrossSectionRMS[iPtRaw],"#sqrt{RMS^{2}+shift^{2}}","f");
    }

    if(!fRelativeVariation) cOutPut[iPtRaw]->cd(1)->DrawFrame(-1.,0.,nFiles,maxrawyield[iPtRaw]*1.5,";Cut set; Raw yield");
    else cOutPut[iPtRaw]->cd(1)->DrawFrame(-1.,0.,nFiles,2.5,";Cut set; Raw yield / Raw yield (central)");
    gRawYieldVsCutSet[iPtRaw]->Draw("PZ");
    if(!fRelativeVariation) cOutPut[iPtRaw]->cd(2)->DrawFrame(-1.,0.,nFiles,maxefficiency[iPtRaw]*1.5,";Cut set; Efficiency");
    else cOutPut[iPtRaw]->cd(2)->DrawFrame(-1.,0.,nFiles,2.5,";Cut set; Efficiency / Efficiency (central)");
    gEffFDVsCutSet[iPtRaw]->Draw("PZ");
    gEffPromptVsCutSet[iPtRaw]->Draw("PZ");
    legEff->Draw("same");
    cOutPut[iPtRaw]->cd(3)->DrawFrame(-1.,0.,nFiles,maxsignif[iPtRaw]*1.5,";Cut set; Significance");
    gSignificanceVsCutSet[iPtRaw]->Draw("PZ");
    cOutPut[iPtRaw]->cd(4)->DrawFrame(-1.,0.,nFiles,maxSoverB[iPtRaw]*1.5,";Cut set; S/B (3#sigma)");
    gSoverBVsCutSet[iPtRaw]->Draw("PZ");
    cOutPut[iPtRaw]->cd(5)->DrawFrame(-1.,maxCross[iPtRaw]/4,nFiles,maxCross[iPtRaw]*1.5,";Cut set; d#it{N}/d#it{p}_{T} (GeV^{-1} #it{c})");
    lCrossSectionCent[iPtRaw]->Draw("same");
    gCrossSectionCent[iPtRaw]->Draw("2");
    gCrossSectionVsCutSet[iPtRaw]->Draw("PZ");
    legCross->Draw("same");
    cOutPut[iPtRaw]->cd(6)->DrawFrame(0.55,0.,1.45,hCrossSectionRatioDist[iPtRaw]->GetMaximum()*1.5,";(d#it{N}/d#it{p}_{T}) / (d#it{N}/d#it{p}_{T})_{central};Entries");
    hCrossSectionRatioDist[iPtRaw]->Draw("same");
    bCrossSectionSyst[iPtRaw]->Draw("same");
    bCrossSectionRMS[iPtRaw]->Draw("same");
    legSyst->Draw("same");
  }
    
  TFile outfile(outfilename.Data(),"recreate");
  for( int iPtRaw=0; iPtRaw<nPtBins; iPtRaw++) {
    cOutPut[iPtRaw]->Write();
    gRawYieldVsCutSet[iPtRaw]->Write();
    gEffFDVsCutSet[iPtRaw]->Write();
    gEffPromptVsCutSet[iPtRaw]->Write();
    gSignificanceVsCutSet[iPtRaw]->Write();
    gSoverBVsCutSet[iPtRaw]->Write();
    gCrossSectionCent[iPtRaw]->Write();
    gCrossSectionVsCutSet[iPtRaw]->Write();
    hCrossSectionRatioDist[iPtRaw]->Write();
  }
  outfile.Close();

  TString outfilenamepdf = outfilename;
  for( int iPtRaw=0; iPtRaw<nPtBins; iPtRaw++) {
    if(iPtRaw==0) outfilenamepdf.ReplaceAll(".root",Form("_ptbin%d.pdf",iPtRaw));
    else outfilenamepdf.ReplaceAll(Form("_ptbin%d.pdf",iPtRaw-1),Form("_ptbin%d.pdf",iPtRaw));
    cOutPut[iPtRaw]->SaveAs(outfilenamepdf.Data());
  }

  return 0;
}

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
