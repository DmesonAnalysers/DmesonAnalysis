//___________________________________________________________________________________//
// Macro for the evaluation of cut-variation systematic uncertainty                  //
// Main Function: PlotCutVariationsOnePtBin                                          //
//___________________________________________________________________________________//

#if !defined (__CINT__) || defined (__CLING__)

#include <iostream>
#include <string>
#include <regex>

#include "yaml-cpp/yaml.h"

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

using namespace std;

//________________________________________________________________________________________________________________
// function prototypes
void PlotCutVariationsOnePtBin(TString cfgFileName = "cfgFile.yml");
void SetStyle();


//________________________________________________________________________________________________________________
// function implementations
void PlotCutVariationsOnePtBin(TString cfgFileName) {
  
    // load inputs
    YAML::Node config = YAML::LoadFile(cfgFileName.Data());
    string inDirName = config["inputs"]["directory"].as<string>();
    string inCommonFileNameRawY = config["inputs"]["commonfilenames"]["rawyield"].as<string>();
    string inCommonFileNameEff = config["inputs"]["commonfilenames"]["efficiency"].as<string>();
    string inCommonFileNameCrossSec = config["inputs"]["commonfilenames"]["crosssection"].as<string>();
    vector<string> cutSetSuffix = config["inputs"]["cutsets"].as<vector<string>>();
    unsigned int nFiles = cutSetSuffix.size();

    string outFileName = config["outfilename"].as<string>();

    double maxChi2 = config["quality"]["maxchisquare"].as<double>();
    double minSignif = config["quality"]["minsignif"].as<double>();
    double minRelSignif = config["quality"]["minrelsignif"].as<double>();
    double minRelEff = config["quality"]["minreleff"].as<double>();
    double maxRelEff = config["quality"]["maxreleff"].as<double>();
    double fillThrRelEff = config["quality"]["fillthrreleff"].as<double>();
    bool fRelativeVariation = static_cast<bool>(config["plots"]["plotrelativevar"].as<int>());
    vector<double> relAssignedSyst = config["plots"]["relassignedsyst"].as<vector<double>>();

    // setting drawing style
    SetStyle();

    TH1D* hCrossSection[nFiles];
    TH1D* hRawYield[nFiles];
    TH1D* hSignificance[nFiles];
    TH1D* hSoverB[nFiles];
    TH1D* hchi2[nFiles];
    TH1D* hEffPrompt[nFiles];
    TH1D* hEffFD[nFiles];

    TString crossSectionTitle = "d#it{N}/d#it{p}_{T}";

    for(unsigned int iFile=0; iFile<nFiles; iFile++) {
        
        TFile* infile_crossec = TFile::Open(Form("%s/%s%s.root",inDirName.data(),Form("crosssec/%s",inCommonFileNameCrossSec.data()),cutSetSuffix[iFile].data()));
        if(!infile_crossec)
            return;
        hCrossSection[iFile] = (TH1D*)infile_crossec->Get("hAAC");
        if(!hCrossSection[iFile]) {
            hCrossSection[iFile] = (TH1D*)infile_crossec->Get("histoSigmaCorr");
            crossSectionTitle = "d#sigma/d#it{p}_{T}";
        }
        if(!hCrossSection[iFile]) {
            hCrossSection[iFile] = (TH1D*)infile_crossec->Get("hCrossSection");
        }
        if(!hCrossSection[iFile]) {
            cerr << "ERROR: cross section histogram not found! Please check it." << endl;
            return;
        }

        hCrossSection[iFile]->SetDirectory(0);
        
        TFile* infile_rawyield = TFile::Open(Form("%s/%s%s.root",inDirName.data(),Form("rawyields/%s",inCommonFileNameRawY.data()),cutSetSuffix[iFile].data()));
        if(!infile_rawyield)
            return;
        hRawYield[iFile] = (TH1D*)infile_rawyield->Get("hRawYields");
        hSignificance[iFile] = (TH1D*)infile_rawyield->Get("hRawYieldsSignificance");
        hSoverB[iFile] = (TH1D*)infile_rawyield->Get("hRawYieldsSoverB");
        hchi2[iFile] = (TH1D*)infile_rawyield->Get("hRawYieldsChiSquare");
        hRawYield[iFile]->SetDirectory(0);
        hSignificance[iFile]->SetDirectory(0);
        hSoverB[iFile]->SetDirectory(0);
        hchi2[iFile]->SetDirectory(0);
        infile_rawyield->Close();

        TFile* infile_eff = TFile::Open(Form("%s/%s%s.root",inDirName.data(),Form("efficiency/%s",inCommonFileNameEff.data()),cutSetSuffix[iFile].data()));
        if(!infile_eff)
            return;
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

    for(int iPt=0; iPt<nPtBins; iPt++) {
        cOutPut[iPt] = new TCanvas(Form("cOutPut_ptbin%d",iPt),"",1920,1080);
        cOutPut[iPt]->Divide(3,2);

        gCrossSectionCent[iPt] = new TGraphErrors(0);
        gCrossSectionCent[iPt]->SetTitle(Form(";Cut set; %s (GeV^{-1} #it{c})", crossSectionTitle.Data()));
        gCrossSectionCent[iPt]->SetName(Form("gCrossSectionCent_ptbin%d",iPt));
        gCrossSectionCent[iPt]->SetLineColor(kRed+1);
        gCrossSectionCent[iPt]->SetFillColorAlpha(kRed+1,0.2);
        gCrossSectionCent[iPt]->SetLineWidth(2);
        
        gCrossSectionVsCutSet[iPt] = new TGraphErrors(0);
        gCrossSectionVsCutSet[iPt]->SetTitle(Form(";Cut set; %s (GeV^{-1} #it{c})", crossSectionTitle.Data()));
        gCrossSectionVsCutSet[iPt]->SetName(Form("gCrossSectionVsCutSet_ptbin%d",iPt));
        gCrossSectionVsCutSet[iPt]->SetMarkerSize(1.);
        gCrossSectionVsCutSet[iPt]->SetMarkerStyle(kFullCircle);
        gCrossSectionVsCutSet[iPt]->SetMarkerColor(kBlack);
        gCrossSectionVsCutSet[iPt]->SetLineColor(kBlack);
        gCrossSectionVsCutSet[iPt]->SetLineWidth(2);

        gRawYieldVsCutSet[iPt] = new TGraphErrors(0);
        gRawYieldVsCutSet[iPt]->SetTitle(";Cut set; raw yield");
        gRawYieldVsCutSet[iPt]->SetName(Form("gRawYieldVsCutSet_ptbin%d",iPt));
        gRawYieldVsCutSet[iPt]->SetMarkerSize(1.);
        gRawYieldVsCutSet[iPt]->SetMarkerStyle(kFullCircle);
        gRawYieldVsCutSet[iPt]->SetMarkerColor(kBlack);
        gRawYieldVsCutSet[iPt]->SetLineColor(kBlack);
        gRawYieldVsCutSet[iPt]->SetLineWidth(2);

        gSignificanceVsCutSet[iPt] = new TGraphErrors(0);
        gSignificanceVsCutSet[iPt]->SetTitle(";Cut set; significance");
        gSignificanceVsCutSet[iPt]->SetName(Form("gSignificanceVsCutSet_ptbin%d",iPt));
        gSignificanceVsCutSet[iPt]->SetMarkerSize(1.);
        gSignificanceVsCutSet[iPt]->SetMarkerStyle(kFullCircle);
        gSignificanceVsCutSet[iPt]->SetMarkerColor(kBlack);
        gSignificanceVsCutSet[iPt]->SetLineColor(kBlack);
        gSignificanceVsCutSet[iPt]->SetLineWidth(2);

        gSoverBVsCutSet[iPt] = new TGraphErrors(0);
        gSoverBVsCutSet[iPt]->SetTitle(";Cut set; S/B (3#sigma)");
        gSoverBVsCutSet[iPt]->SetName(Form("gSoverBVsCutSet_ptbin%d",iPt));
        gSoverBVsCutSet[iPt]->SetMarkerSize(1.);
        gSoverBVsCutSet[iPt]->SetMarkerStyle(kFullCircle);
        gSoverBVsCutSet[iPt]->SetMarkerColor(kBlack);
        gSoverBVsCutSet[iPt]->SetLineColor(kBlack);
        gSoverBVsCutSet[iPt]->SetLineWidth(2);

        gEffPromptVsCutSet[iPt] = new TGraphErrors(0);
        gEffPromptVsCutSet[iPt]->SetTitle(";Cut set; Efficiency");
        gEffPromptVsCutSet[iPt]->SetName(Form("gEffPromptVsCutSet_ptbin%d",iPt));
        gEffPromptVsCutSet[iPt]->SetMarkerSize(1.);
        gEffPromptVsCutSet[iPt]->SetMarkerStyle(kFullCircle);
        gEffPromptVsCutSet[iPt]->SetMarkerColor(kRed);
        gEffPromptVsCutSet[iPt]->SetLineColor(kRed);
        gEffPromptVsCutSet[iPt]->SetLineWidth(2);

        gEffFDVsCutSet[iPt] = new TGraphErrors(0);
        gEffFDVsCutSet[iPt]->SetTitle(";Cut set; Efficiency");
        gEffFDVsCutSet[iPt]->SetName(Form("gEffFDVsCutSet_ptbin%d",iPt));
        gEffFDVsCutSet[iPt]->SetMarkerSize(1.);
        gEffFDVsCutSet[iPt]->SetMarkerStyle(kFullSquare);
        gEffFDVsCutSet[iPt]->SetMarkerColor(kBlue);
        gEffFDVsCutSet[iPt]->SetLineColor(kBlue);
        gEffFDVsCutSet[iPt]->SetLineWidth(2);

        hCrossSectionRatioDist[iPt] = new TH1D(Form("hCrossSectionRatioDist_ptbin%d",iPt), Form(";(%s) / (%s)_{central};Entries", crossSectionTitle.Data(), crossSectionTitle.Data()),60,0.45,1.55);
        hCrossSectionRatioDist[iPt]->SetLineColor(kBlack);
        hCrossSectionRatioDist[iPt]->SetLineWidth(2);

        maxrawyield[iPt]=-1;
        maxefficiency[iPt]=-1;
        maxsignif[iPt]=-1;
        maxSoverB[iPt]=-1;
        maxCross[iPt]=-1;

        for(unsigned int iFile=0; iFile<nFiles; iFile++) {
            if(maxrawyield[iPt]<hRawYield[iFile]->GetBinContent(iPt+1))
                maxrawyield[iPt] = hRawYield[iFile]->GetBinContent(iPt+1);
            if(maxefficiency[iPt]<hEffFD[iFile]->GetBinContent(iPt+1))
                maxefficiency[iPt] = hEffFD[iFile]->GetBinContent(iPt+1);
            if(maxsignif[iPt]<hSignificance[iFile]->GetBinContent(iPt+1))
                maxsignif[iPt] = hSignificance[iFile]->GetBinContent(iPt+1);
            if(maxSoverB[iPt]<hSoverB[iFile]->GetBinContent(iPt+1))
                maxSoverB[iPt] = hSoverB[iFile]->GetBinContent(iPt+1);
            if(maxCross[iPt]<hCrossSection[iFile]->GetBinContent(iPt+1))
                maxCross[iPt] = hCrossSection[iFile]->GetBinContent(iPt+1);
            
            int crossptbin = hCrossSection[iFile]->GetXaxis()->FindBin(hRawYield[iFile]->GetBinCenter(iPt+1));
            int effptbin = hEffPrompt[iFile]->GetXaxis()->FindBin(hRawYield[iFile]->GetBinCenter(iPt+1));

            if(iFile==0) {
                lCrossSectionCent[iPt] = new TLine(0.,hCrossSection[iFile]->GetBinContent(crossptbin),nFiles-1,hCrossSection[iFile]->GetBinContent(crossptbin));
                lCrossSectionCent[iPt]->SetLineWidth(2);
                lCrossSectionCent[iPt]->SetLineColor(kRed+1);

                gCrossSectionCent[iPt]->SetPoint(iFile,nFiles/2,hCrossSection[iFile]->GetBinContent(crossptbin));
                gCrossSectionCent[iPt]->SetPointError(iFile,nFiles/2,hCrossSection[iFile]->GetBinContent(crossptbin)*relAssignedSyst[iPt]);
            }
          
            gCrossSectionVsCutSet[iPt]->SetPoint(iFile,iFile,hCrossSection[iFile]->GetBinContent(crossptbin));
            gCrossSectionVsCutSet[iPt]->SetPointError(iFile,0.5,hCrossSection[iFile]->GetBinError(crossptbin));

            if(!fRelativeVariation) {
                gRawYieldVsCutSet[iPt]->SetPoint(iFile,iFile,hRawYield[iFile]->GetBinContent(iPt+1));
                gRawYieldVsCutSet[iPt]->SetPointError(iFile,0.5,hRawYield[iFile]->GetBinError(iPt+1));

                gEffPromptVsCutSet[iPt]->SetPoint(iFile,iFile,hEffPrompt[iFile]->GetBinContent(effptbin));
                gEffPromptVsCutSet[iPt]->SetPointError(iFile,0.5,hEffPrompt[iFile]->GetBinError(effptbin));

                gEffFDVsCutSet[iPt]->SetPoint(iFile,iFile,hEffFD[iFile]->GetBinContent(effptbin));
                gEffFDVsCutSet[iPt]->SetPointError(iFile,0.5,hEffFD[iFile]->GetBinError(effptbin));
            }
            else {
                gRawYieldVsCutSet[iPt]->SetPoint(iFile,iFile,hRawYield[iFile]->GetBinContent(iPt+1)/hRawYield[0]->GetBinContent(iPt+1));
                gRawYieldVsCutSet[iPt]->SetPointError(iFile,0.5,0);

                gEffPromptVsCutSet[iPt]->SetPoint(iFile,iFile,hEffPrompt[iFile]->GetBinContent(effptbin)/hEffPrompt[0]->GetBinContent(effptbin));
                gEffPromptVsCutSet[iPt]->SetPointError(iFile,0.5,0);

                gEffFDVsCutSet[iPt]->SetPoint(iFile,iFile,hEffFD[iFile]->GetBinContent(effptbin)/hEffFD[0]->GetBinContent(effptbin));
                gEffFDVsCutSet[iPt]->SetPointError(iFile,0.5,0);
            }

            gSignificanceVsCutSet[iPt]->SetPoint(iFile,iFile,hSignificance[iFile]->GetBinContent(iPt+1));
            gSignificanceVsCutSet[iPt]->SetPointError(iFile,0.5,hSignificance[iFile]->GetBinError(iPt+1));

            gSoverBVsCutSet[iPt]->SetPoint(iFile,iFile,hSoverB[iFile]->GetBinContent(iPt+1));
            gSoverBVsCutSet[iPt]->SetPointError(iFile,0.5,hSoverB[iFile]->GetBinError(iPt+1));      
    
            if(iFile!=0) {
                if(hchi2[iFile]->GetBinContent(iPt+1) > maxChi2)
                    continue;
                if(TMath::Abs(1-hEffPrompt[iFile]->GetBinContent(effptbin)/hEffPrompt[0]->GetBinContent(effptbin)) < fillThrRelEff && 
                   TMath::Abs(1-hEffFD[iFile]->GetBinContent(iPt+1)/hEffFD[0]->GetBinContent(iPt+1)) < fillThrRelEff)
                    continue;
                if((hEffPrompt[iFile]->GetBinContent(effptbin)/hEffPrompt[0]->GetBinContent(effptbin)) < minRelEff)
                    continue;
                if((hEffPrompt[iFile]->GetBinContent(effptbin)/hEffPrompt[0]->GetBinContent(effptbin)) > maxRelEff)
                    continue;
                if(hSignificance[iFile]->GetBinContent(iPt+1) < minSignif)
                    continue;
                if((hSignificance[iFile]->GetBinContent(iPt+1)/hSignificance[0]->GetBinContent(iPt+1)) < minRelSignif)
                    continue;
                hCrossSectionRatioDist[iPt]->Fill(hCrossSection[iFile]->GetBinContent(crossptbin)/hCrossSection[0]->GetBinContent(crossptbin));
            }
        }

        bCrossSectionSyst[iPt] = new TBox(1-relAssignedSyst[iPt],0.,1+relAssignedSyst[iPt],hCrossSectionRatioDist[iPt]->GetMaximum());
        bCrossSectionSyst[iPt]->SetLineColor(kRed+1);
        bCrossSectionSyst[iPt]->SetFillColorAlpha(kRed+1,0.2);

        bCrossSectionRMS[iPt] = new TBox(1-TMath::Sqrt(hCrossSectionRatioDist[iPt]->GetRMS()*hCrossSectionRatioDist[iPt]->GetRMS()+(1-hCrossSectionRatioDist[iPt]->GetMean())*(1-hCrossSectionRatioDist[iPt]->GetMean())),0.,1+TMath::Sqrt(hCrossSectionRatioDist[iPt]->GetRMS()*hCrossSectionRatioDist[iPt]->GetRMS()+(1-hCrossSectionRatioDist[iPt]->GetMean())*(1-hCrossSectionRatioDist[iPt]->GetMean())),hCrossSectionRatioDist[iPt]->GetMaximum());
        bCrossSectionRMS[iPt]->SetLineColor(kBlue+1);
        bCrossSectionRMS[iPt]->SetFillColorAlpha(kBlue+1,0.2);

        if(iPt==0) {
            legEff->AddEntry(gEffPromptVsCutSet[iPt],"Prompt","lpe");
            legEff->AddEntry(gEffFDVsCutSet[iPt],"Feed-down","lpe");
            legCross->AddEntry(gCrossSectionCent[iPt],"Central value #pm assigned syst","fl");
            legSyst->AddEntry(bCrossSectionSyst[iPt],"Assigned syst","f");
            legSyst->AddEntry(bCrossSectionRMS[iPt],"#sqrt{RMS^{2}+shift^{2}}","f");
        }

        if(!fRelativeVariation) 
            cOutPut[iPt]->cd(1)->DrawFrame(-1.,0.,nFiles,maxrawyield[iPt]*1.5,";Cut set; Raw yield");
        else 
            cOutPut[iPt]->cd(1)->DrawFrame(-1.,0.,nFiles,2.5,";Cut set; Raw yield / Raw yield (central)");
        gRawYieldVsCutSet[iPt]->Draw("PZ");
        if(!fRelativeVariation) 
            cOutPut[iPt]->cd(2)->DrawFrame(-1.,0.,nFiles,maxefficiency[iPt]*1.5,";Cut set; Efficiency");
        else 
            cOutPut[iPt]->cd(2)->DrawFrame(-1.,0.,nFiles,2.5,";Cut set; Efficiency / Efficiency (central)");
        gEffFDVsCutSet[iPt]->Draw("PZ");
        gEffPromptVsCutSet[iPt]->Draw("PZ");
        legEff->Draw("same");
        cOutPut[iPt]->cd(3)->DrawFrame(-1.,0.,nFiles,maxsignif[iPt]*1.5,";Cut set; Significance");
        gSignificanceVsCutSet[iPt]->Draw("PZ");
        cOutPut[iPt]->cd(4)->DrawFrame(-1.,0.,nFiles,maxSoverB[iPt]*1.5,";Cut set; S/B (3#sigma)");
        gSoverBVsCutSet[iPt]->Draw("PZ");
        cOutPut[iPt]->cd(5)->DrawFrame(-1.,maxCross[iPt]/4,nFiles,maxCross[iPt]*1.5, Form(";Cut set; %s (GeV^{-1} #it{c})", crossSectionTitle.Data()));
        lCrossSectionCent[iPt]->Draw("same");
        gCrossSectionCent[iPt]->Draw("2");
        gCrossSectionVsCutSet[iPt]->Draw("PZ");
        legCross->Draw("same");
        cOutPut[iPt]->cd(6)->DrawFrame(0.55,0.,1.45,hCrossSectionRatioDist[iPt]->GetMaximum()*1.5, Form(";(%s) / (%s)_{central};Entries", crossSectionTitle.Data(), crossSectionTitle.Data()));
        hCrossSectionRatioDist[iPt]->Draw("same");
        bCrossSectionSyst[iPt]->Draw("same");
        bCrossSectionRMS[iPt]->Draw("same");
        legSyst->Draw("same");
    }
    
    TFile outfile(outFileName.data(),"recreate");
    for(int iPt=0; iPt<nPtBins; iPt++) {
        cOutPut[iPt]->Write();
        gRawYieldVsCutSet[iPt]->Write();
        gEffFDVsCutSet[iPt]->Write();
        gEffPromptVsCutSet[iPt]->Write();
        gSignificanceVsCutSet[iPt]->Write();
        gSoverBVsCutSet[iPt]->Write();
        gCrossSectionCent[iPt]->Write();
        gCrossSectionVsCutSet[iPt]->Write();
        hCrossSectionRatioDist[iPt]->Write();
    }
    outfile.Close();

    outFileName = regex_replace(outFileName, regex(".root"), ".pdf");
    cOutPut[0]->SaveAs(Form("%s[", outFileName.data()));
    for(int iPt=0; iPt<nPtBins; iPt++)
        cOutPut[iPt]->SaveAs(outFileName.data());
    cOutPut[nPtBins-1]->SaveAs(Form("%s]", outFileName.data()));
}


//________________________________________________________________________________________________________________
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
    gStyle->SetTitleSize(0.05,"xy");
    gStyle->SetLabelSize(0.05,"xy");
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadRightMargin(0.12);
    TGaxis::SetMaxDigits(3);
}
