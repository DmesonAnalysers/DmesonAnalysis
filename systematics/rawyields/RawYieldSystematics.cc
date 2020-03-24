//___________________________________________________________________________________//
// Macro for the evaluation of raw-yield extraction systematic uncertainty           //
// Main Function: RawYieldSystematics                                                //
//___________________________________________________________________________________//

#if !defined(__CINT__) || defined(__CLING__)

#include <iostream>
#include <string>
#include <regex>

#include "yaml-cpp/yaml.h"

#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>

#include "AliHFInvMassFitter.h"
#include "AliHFInvMassMultiTrialFit.h"
#include "AliVertexingHFUtils.h"

#endif

//________________________________________________________________________________________________________________
// function prototypes
int RawYieldSystematics(TString cfgFileName = "cfgFile.yml");
double min(TH1F *histo);
double max(TH1F *histo);
int LoadRefFiles(string refFileName, string refFileNameMC,
                 TH1F *&hRawYieldRef, TH1F *&hSigmaRef,
                 TH1F *&hMeanRef, TH1F *&hSigmaMC);
void SetStyle();


//________________________________________________________________________________________________________________
// function implementations
int RawYieldSystematics(TString cfgFileName) {

    // load inputs
    YAML::Node config = YAML::LoadFile(cfgFileName.Data());
    string refFileName = config["reffilenames"]["data"].as<string>();
    string refFileNameMC = config["reffilenames"]["MC"].as<string>();
    string outFileName = config["outfilename"].as<string>();
    vector<double> mins = config["multitrial"]["mins"].as<vector<double>>();
    vector<double> maxs = config["multitrial"]["maxs"].as<vector<double>>();
    vector<int> rebins = config["multitrial"]["rebins"].as<vector<int>>();
    vector<string> signalFuncs = config["multitrial"]["sgnfuncs"].as<vector<string>>();
    vector<string> bkgFuncs = config["multitrial"]["bkgfuncs"].as<vector<string>>();
    vector<string> sigmaOpt = config["multitrial"]["sigma"].as<vector<string>>();
    const unsigned int nSigmaConf = sigmaOpt.size();
    vector<string> meanOpt = config["multitrial"]["mean"].as<vector<string>>();
    const unsigned int nMeanConf = meanOpt.size();
    vector<double> nSigmaBinCounting = config["multitrial"]["bincounting"]["nsigma"].as<vector<double>>();
    const unsigned int nBinCounting = nSigmaBinCounting.size();
    double minChi2 = config["quality"]["chisquare"]["min"].as<double>();
    double maxChi2 = config["quality"]["chisquare"]["max"].as<double>();

    string mesonName = config["meson"].as<string>();    
    double massDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
    double massDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    double mass = -1;
    if(mesonName == "Ds")
        mass = massDs;
    else if(mesonName == "Dplus")
        mass = massDplus;
    else {
        std::cerr << "ERROR: you must specify if it is D+ or Ds+! Exit." << std::endl;
        return -1;
    }

    TH1F *hRawYieldRef = nullptr;
    TH1F *hSigmaRef = nullptr;
    TH1F *hMeanRef = nullptr;
    TH1F *hSigmaMC = nullptr;
    int loadref = LoadRefFiles(refFileName, refFileNameMC, hRawYieldRef, hSigmaRef, hMeanRef, hSigmaMC);
    if (loadref > 0) {
        std::cerr << "ERROR: missing information in reference files! Check them please." << std::endl;
        return loadref;
    }

    const int nPtBins = hRawYieldRef->GetNbinsX();
    const int nPtLims = nPtBins + 1;
    double PtLims[nPtLims];
    TH1D *hMass[nPtBins];
    TFile *infile = TFile::Open(refFileName.data());
    for (int iPt = 0; iPt < nPtBins; iPt++)
        PtLims[iPt] = hRawYieldRef->GetBinLowEdge(iPt+1);
    PtLims[nPtBins] = hRawYieldRef->GetBinLowEdge(nPtBins) + hRawYieldRef->GetBinWidth(nPtBins);
    for (int iPt = 0; iPt < nPtBins; iPt++)
        hMass[iPt] = (TH1D*)infile->Get(Form("hMass_%0.f_%0.f", PtLims[iPt] * 10, PtLims[iPt+1] * 10));
    
    // outputs
    TNtuple* ntupleFit[nPtBins];
    TNtuple* ntupleBinC[nPtBins];
    TH1F *hRawYield[nPtBins];
    TH1F *hBinCount[nBinCounting][nPtBins];
    TH1F *hSigma[nPtBins];
    TH1F *hMean[nPtBins];
    TH1F *hChiSquare[nPtBins];
    TH1F *hRawYieldVsTrial[nPtBins];
    TH1F *hBinCountVsTrial[nBinCounting][nPtBins];
    TH1F *hSigmaVsTrial[nPtBins];
    TH1F *hMeanVsTrial[nPtBins];
    TH1F *hChiSquareVsTrial[nPtBins];

    TLine *lRawRef[nPtBins];
    TLine *lSigmaRef[nPtBins];
    TLine *lMeanRef[nPtBins];
    TLine *lRawRefVsTrial[nPtBins];
    TLine *lSigmaRefVsTrial[nPtBins];
    TLine *lMeanRefVsTrial[nPtBins];

    TPaveText *stats[nPtBins];
    TPaveText *statsbc[nBinCounting][nPtBins];
    TLegend *leg = new TLegend(0.15, 0.65, 0.35, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.05);

    // setting drawing style
    SetStyle();
    int binCountColors[] = {kOrange + 1, kGreen + 2, kRed + 1};
    int binCountColorsMarker[] = {kBlack, kCyan, kBlack};

    for (int iPt = 0; iPt < nPtBins; iPt++) {

        // configure multi-trial
        AliHFInvMassMultiTrialFit multiTrial;
        multiTrial.SetUseLogLikelihoodFit();
        multiTrial.SetMass(mass);
        if(mesonName == "Ds")
            multiTrial.IncludeSecondGausPeak(massDplus, true, 0.9*hSigmaMC->GetBinContent(iPt+1), true);

        multiTrial.ConfigureLowLimFitSteps(mins.size(), mins.data());
        multiTrial.ConfigureUpLimFitSteps(maxs.size(), maxs.data());
        multiTrial.ConfigureRebinSteps(rebins.size(), rebins.data());
        multiTrial.ConfigurenSigmaBinCSteps(nSigmaBinCounting.size(), nSigmaBinCounting.data());
        
        multiTrial.SetUseExpoBackground(false);
        multiTrial.SetUseLinBackground(false);
        multiTrial.SetUsePol2Background(false);
        multiTrial.SetUsePol3Background(false);
        multiTrial.SetUsePol4Background(false);
        multiTrial.SetUsePol5Background(false);
        multiTrial.SetUsePowerLawBackground(false);
        multiTrial.SetUsePowerLawTimesExpoBackground(false);
        for(auto func: bkgFuncs) {
            if(func == "kExpo")
                multiTrial.SetUseExpoBackground();
            else if(func == "kLin")
                multiTrial.SetUseLinBackground();
            else if(func == "kPol2")
                multiTrial.SetUsePol2Background();
            else if(func == "kPol3")
                multiTrial.SetUsePol3Background();
            else if(func == "kPol4")
                multiTrial.SetUsePol4Background();
            else if(func == "kPol5")
                multiTrial.SetUsePol5Background();
        }

        for(auto func: signalFuncs) {
            if(func == "k2Gaus")
                multiTrial.SetUse2GausSignal();
        }

        multiTrial.SetUseFixedMeanFreeS(false);
        for(auto opt: meanOpt) {
            if(opt == "kFixed")
                multiTrial.SetUseFixedMeanFreeS();
        }

        multiTrial.SetUseFixSigUpFreeMean(false);
        multiTrial.SetUseFixSigDownFreeMean(false);
        multiTrial.SetUseFreeS(false);
        multiTrial.SetUseFixSigFreeMean(false);
        multiTrial.SetUseFixSigFixMean(false);
        multiTrial.SetUseFixSigFixMeanUp(false);
        multiTrial.SetUseFixSigFixMeanDown(false);
        multiTrial.SetUseFreeSigFixMeanUp(false);
        multiTrial.SetUseFreeSigFixMeanDown(false);
        multiTrial.SetUseFixSigVarWithFixMean(false);
        for(auto opt: sigmaOpt) {
            if(opt == "kFree")
                multiTrial.SetUseFreeS();
            else if(opt == "kFixed") {
                multiTrial.SetSigmaGaussMC(hSigmaMC->GetBinContent(iPt+1));
                multiTrial.SetUseFixSigFreeMean();
            }
            else if(opt == "kFixedMinus10Perc") {
                multiTrial.SetSigmaMCVariation(0.10);
                multiTrial.SetUseFixSigDownFreeMean();
            }
            else if(opt == "kFixedPlus10Perc") {
                multiTrial.SetSigmaMCVariation(0.10);
                multiTrial.SetUseFixSigUpFreeMean();
            }
            else if(opt == "kFixedMinus15Perc") {
                multiTrial.SetSigmaMCVariation(0.15);
                multiTrial.SetUseFixSigDownFreeMean();
            }
            else if(opt == "kFixedPlus15Perc") {
                multiTrial.SetSigmaMCVariation(0.15);
                multiTrial.SetUseFixSigUpFreeMean();
            }
            else if(opt == "kFixedMinus20Perc") {
                multiTrial.SetSigmaMCVariation(0.15);
                multiTrial.SetUseFixSigDownFreeMean();
            }
            else if(opt == "kFixedPlus20Perc") {
                multiTrial.SetSigmaMCVariation(0.15);
                multiTrial.SetUseFixSigUpFreeMean();
            }
            else if(opt == "kFixedMinusUnc") {
                multiTrial.SetSigmaMCVariation(hSigmaMC->GetBinError(iPt+1) / hSigmaMC->GetBinContent(iPt+1));
                multiTrial.SetUseFixSigDownFreeMean();
            }
            else if(opt == "kFixedPlusUnc") {
                multiTrial.SetSigmaMCVariation(hSigmaMC->GetBinError(iPt+1) / hSigmaMC->GetBinContent(iPt+1));
                multiTrial.SetUseFixSigUpFreeMean();
            }
        }

        // perform multi-trial
        std::cout << "\n\n*****************************************************" << std::endl;
        std::cout << Form("Perform multi-trial for %0.f < pT < %0.f GeV/c", PtLims[iPt], PtLims[iPt+1]) << std::endl;
        multiTrial.DoMultiTrials(hMass[iPt]);
        ntupleFit[iPt] = (TNtuple *)(multiTrial.GetNtupleMultiTrials())->Clone(Form("ntupleFit_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]));
        ntupleBinC[iPt] = (TNtuple *)(multiTrial.GetNtupleBinCounting())->Clone(Form("ntupleBinC_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]));
        ntupleFit[iPt]->SetDirectory(0);
        ntupleBinC[iPt]->SetDirectory(0);

        // fill outputs
        unsigned int minSigmaConf = ntupleFit[iPt]->GetMinimum("confsig");
        unsigned int maxSigmaConf = ntupleFit[iPt]->GetMaximum("confsig");
        unsigned int minMeanConf = ntupleFit[iPt]->GetMinimum("confmean");
        unsigned int maxMeanConf = ntupleFit[iPt]->GetMaximum("confmean");

        int nTrialsFit = ntupleFit[iPt]->GetEntries(Form("chi2 >= %f && chi2 <= %f", minChi2, maxChi2));
        int nTrialsBinC = ntupleBinC[iPt]->GetEntries(Form("chi2 >= %f && chi2 <= %f && confsig == %d && confmean == %d",
                                                           minChi2, maxChi2, minSigmaConf, minMeanConf)); // keep only one configuration of sigma and mean for BC
        int nTrials = nTrialsFit + nTrialsBinC;

        double rawMin = ntupleFit[iPt]->GetMinimum("rawy") * 0.5;
        double rawMax = ntupleFit[iPt]->GetMaximum("rawy") * 1.5;
        if (loadref != 1 && loadref != 2 && loadref != 4) {
            lRawRefVsTrial[iPt] = new TLine(-0.5, hRawYieldRef->GetBinContent(iPt+1), nTrials - 0.5, hRawYieldRef->GetBinContent(iPt+1));
            lRawRefVsTrial[iPt]->SetLineColor(kRed);
            lRawRefVsTrial[iPt]->SetLineWidth(2);
        }

        double sigmaMin = ntupleFit[iPt]->GetMinimum("sigma") * 0.5;
        double sigmaMax = ntupleFit[iPt]->GetMaximum("sigma") * 1.5;
        if (loadref != 1 && loadref != 2 && loadref != 5) {
            lSigmaRefVsTrial[iPt] = new TLine(-0.5, hSigmaRef->GetBinContent(iPt+1), nTrials - 0.5, hSigmaRef->GetBinContent(iPt+1));
            lSigmaRefVsTrial[iPt]->SetLineColor(kRed);
            lSigmaRefVsTrial[iPt]->SetLineWidth(2);
        }

        double meanMin = ntupleFit[iPt]->GetMinimum("mean") * (1 - 0.005);
        double meanMax = ntupleFit[iPt]->GetMaximum("mean") * (1 + 0.005);
        if (loadref != 1 && loadref != 2 && loadref != 6) {
            lMeanRefVsTrial[iPt] = new TLine(-0.5, hMeanRef->GetBinContent(iPt+1), nTrials - 0.5, hMeanRef->GetBinContent(iPt+1));
            lMeanRefVsTrial[iPt]->SetLineColor(kRed);
            lMeanRefVsTrial[iPt]->SetLineWidth(2);
        }

        int nBins = 100;

        hRawYield[iPt] = new TH1F(Form("hRawYield_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]), 
                                  Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};raw yield;entries", PtLims[iPt], PtLims[iPt+1]), nBins, rawMin, rawMax);
        hRawYield[iPt]->SetFillStyle(3004);
        hRawYield[iPt]->SetLineWidth(2);
        hRawYield[iPt]->SetLineColor(kBlue + 1);
        hRawYield[iPt]->SetFillColor(kBlue + 1);
        ntupleFit[iPt]->Draw(Form("rawy>>hRawYield_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]), Form("chi2 > %f && chi2 < %f", minChi2, maxChi2));

        hSigma[iPt] = new TH1F(Form("hSigma_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]), 
                               Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};width (GeV/c^{2});Entries", PtLims[iPt], PtLims[iPt+1]), nBins, sigmaMin, sigmaMax);
        hSigma[iPt]->SetFillStyle(3004);
        hSigma[iPt]->SetLineWidth(2);
        hSigma[iPt]->SetLineColor(kBlue + 1);
        hSigma[iPt]->SetFillColor(kBlue + 1);
        ntupleFit[iPt]->Draw(Form("sigma>>hSigma_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]), Form("chi2 > %f && chi2 < %f", minChi2, maxChi2));

        hMean[iPt] = new TH1F(Form("hMean_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]),
                              Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};mean (GeV/c^{2});entries", PtLims[iPt], PtLims[iPt+1]), nBins, meanMin, meanMax);
        hMean[iPt]->SetFillStyle(3004);
        hMean[iPt]->SetLineWidth(2);
        hMean[iPt]->SetLineColor(kBlue + 1);
        hMean[iPt]->SetFillColor(kBlue + 1);
        ntupleFit[iPt]->Draw(Form("mean>>hMean_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]), Form("chi2 > %f && chi2 < %f", minChi2, maxChi2));

        hChiSquare[iPt] = new TH1F(Form("hChiSquare_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]),
                                   Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};#chi^{2}/ndf;entries", PtLims[iPt], PtLims[iPt+1]), nBins, 0., maxChi2);
        hChiSquare[iPt]->SetFillStyle(3004);
        hChiSquare[iPt]->SetLineWidth(2);
        hChiSquare[iPt]->SetLineColor(kBlue + 1);
        hChiSquare[iPt]->SetFillColor(kBlue + 1);
        ntupleFit[iPt]->Draw(Form("chi2>>hChiSquare_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]), Form("chi2 > %f && chi2 < %f", minChi2, maxChi2));

        hRawYieldVsTrial[iPt] = new TH1F(Form("hRawYieldVsTrial_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]),
                                         Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};trial # ;raw yield", PtLims[iPt], PtLims[iPt+1]), nTrials, -0.5, nTrials - 0.5);
        hRawYieldVsTrial[iPt]->GetYaxis()->SetRangeUser(0., rawMax*1.5);
        hRawYieldVsTrial[iPt]->SetMarkerStyle(20);
        hRawYieldVsTrial[iPt]->SetLineWidth(2);
        hRawYieldVsTrial[iPt]->SetMarkerSize(0.5);
        hRawYieldVsTrial[iPt]->SetLineColor(kBlue + 1);
        hRawYieldVsTrial[iPt]->SetMarkerColor(kBlack);

        hSigmaVsTrial[iPt] = new TH1F(Form("hMeanVsTrial_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]),
                                      Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};trial # ;width (GeV/c^{2})", PtLims[iPt], PtLims[iPt+1]), nTrials, -0.5, nTrials - 0.5);
        hSigmaVsTrial[iPt]->GetYaxis()->SetRangeUser(sigmaMin, sigmaMax);
        hSigmaVsTrial[iPt]->SetMarkerStyle(20);
        hSigmaVsTrial[iPt]->SetLineWidth(2);
        hSigmaVsTrial[iPt]->SetMarkerSize(0.5);
        hSigmaVsTrial[iPt]->SetLineColor(kBlue + 1);
        hSigmaVsTrial[iPt]->SetMarkerColor(kBlack);

        hMeanVsTrial[iPt] = new TH1F(Form("hSigmaVsTrial_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]),
                                     Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};trial # ;mean (GeV/c^{2})", PtLims[iPt], PtLims[iPt+1]), nTrials, -0.5, nTrials - 0.5);
        hMeanVsTrial[iPt]->GetYaxis()->SetRangeUser(meanMin, meanMax);
        hMeanVsTrial[iPt]->SetMarkerStyle(20);
        hMeanVsTrial[iPt]->SetLineWidth(2);
        hMeanVsTrial[iPt]->SetMarkerSize(0.5);
        hMeanVsTrial[iPt]->SetLineColor(kBlue + 1);
        hMeanVsTrial[iPt]->SetMarkerColor(kBlack);

        hChiSquareVsTrial[iPt] = new TH1F(Form("hChiSquareVsTrial_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]),
                                          Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};trial # ;#chi^{2}/ndf", PtLims[iPt], PtLims[iPt+1]), nTrials, -0.5, nTrials - 0.5);
        hChiSquareVsTrial[iPt]->GetYaxis()->SetRangeUser(0, maxChi2);
        hChiSquareVsTrial[iPt]->SetMarkerStyle(20);
        hChiSquareVsTrial[iPt]->SetLineWidth(2);
        hChiSquareVsTrial[iPt]->SetMarkerSize(0.5);
        hChiSquareVsTrial[iPt]->SetLineColor(kBlue + 1);
        hChiSquareVsTrial[iPt]->SetMarkerColor(kBlack);

        float rawYield, rawYielUnc, sigma, sigmaUnc, pos, posUnc, chiSquare, confSigma, confMean;
        ntupleFit[iPt]->SetBranchAddress("rawy", &rawYield);
        ntupleFit[iPt]->SetBranchAddress("erawy", &rawYielUnc);
        ntupleFit[iPt]->SetBranchAddress("sigma", &sigma);
        ntupleFit[iPt]->SetBranchAddress("esigma", &sigmaUnc);
        ntupleFit[iPt]->SetBranchAddress("mean", &pos);
        ntupleFit[iPt]->SetBranchAddress("emean", &posUnc);
        ntupleFit[iPt]->SetBranchAddress("chi2", &chiSquare);
        ntupleFit[iPt]->SetBranchAddress("confsig", &confSigma);
        ntupleFit[iPt]->SetBranchAddress("confmean", &confMean);

        int iBin = 1;
        for(unsigned int iSigma = minSigmaConf; iSigma <= minSigmaConf+(maxSigmaConf-minSigmaConf); iSigma++) {
            for(unsigned int iMean = minMeanConf; iMean <= minMeanConf+(maxMeanConf-minMeanConf); iMean++) {
                for(int iEntry = 0; iEntry < ntupleFit[iPt]->GetEntriesFast(); iEntry++) {
                    ntupleFit[iPt]->GetEntry(iEntry);
                    if(chiSquare > maxChi2 || chiSquare < minChi2 || TMath::Abs(confSigma-iSigma) > 0.001 || TMath::Abs(confMean-iMean) > 0.001)
                        continue;
                    hRawYieldVsTrial[iPt]->SetBinContent(iBin, rawYield);
                    hRawYieldVsTrial[iPt]->SetBinError(iBin, rawYielUnc);
                    hSigmaVsTrial[iPt]->SetBinContent(iBin, sigma);
                    hSigmaVsTrial[iPt]->SetBinError(iBin, sigmaUnc);
                    hMeanVsTrial[iPt]->SetBinContent(iBin, pos);
                    hMeanVsTrial[iPt]->SetBinError(iBin, posUnc);
                    hChiSquareVsTrial[iPt]->SetBinContent(iBin, chiSquare);
                    iBin++;
                }
            }
        }
        
        for (unsigned int iBinCount = 0; iBinCount < nBinCounting; iBinCount++) {
            hBinCount[iBinCount][iPt] = new TH1F(Form("hBinCount_%0.1fsigma_pT_%0.f-%0.f", nSigmaBinCounting[iBinCount], PtLims[iPt], PtLims[iPt+1]), 
                                                Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};raw yield;entries", PtLims[iPt], PtLims[iPt+1]), nBins, rawMin, rawMax);
            hBinCount[iBinCount][iPt]->SetFillStyle(3004);
            hBinCount[iBinCount][iPt]->SetLineWidth(2);
            hBinCount[iBinCount][iPt]->SetLineColor(binCountColors[iBinCount]);
            hBinCount[iBinCount][iPt]->SetFillColor(binCountColors[iBinCount]);
            ntupleBinC[iPt]->Draw(Form("rawyBC1>>hBinCount_%0.1fsigma_pT_%0.f-%0.f", nSigmaBinCounting[iBinCount], PtLims[iPt], PtLims[iPt+1]),
                                  Form("chi2 > %f && chi2 < %f && nSigmaBC == %f", minChi2, maxChi2, nSigmaBinCounting[iBinCount]));

            hBinCountVsTrial[iBinCount][iPt] = new TH1F(Form("hBinCountVsTrial_%0.1fsigma_pT_%0.f-%0.f", nSigmaBinCounting[iBinCount], PtLims[iPt], PtLims[iPt+1]),
                                                        Form("%0.f < #it{p}_{T} < %0.f GeV/#it{c};trial # ;raw yield", PtLims[iPt], PtLims[iPt+1]), nTrials, -0.5, nTrials - 0.5);
            hBinCountVsTrial[iBinCount][iPt]->GetYaxis()->SetRangeUser(0., rawMax*1.5);
            hBinCountVsTrial[iBinCount][iPt]->SetMarkerStyle(20);
            hBinCountVsTrial[iBinCount][iPt]->SetLineWidth(2);
            hBinCountVsTrial[iBinCount][iPt]->SetMarkerSize(0.5);
            hBinCountVsTrial[iBinCount][iPt]->SetLineColor(binCountColors[iBinCount]);
            hBinCountVsTrial[iBinCount][iPt]->SetMarkerColor(binCountColorsMarker[iBinCount]);

            float rawYieldBC, rawYielBCUnc, chiSquareBC, nSigmaBC, confSigmaBC, confMeanBC;
            ntupleBinC[iPt]->SetBranchAddress("rawyBC1", &rawYieldBC);
            ntupleBinC[iPt]->SetBranchAddress("erawyBC1", &rawYielBCUnc);
            ntupleBinC[iPt]->SetBranchAddress("nSigmaBC", &nSigmaBC);
            ntupleBinC[iPt]->SetBranchAddress("chi2", &chiSquareBC);
            ntupleBinC[iPt]->SetBranchAddress("confsig", &confSigmaBC);
            ntupleBinC[iPt]->SetBranchAddress("confmean", &confMeanBC);

            int Counter = 0;
            for(int iEntry = 0; iEntry < ntupleBinC[iPt]->GetEntriesFast(); iEntry++) {
                ntupleBinC[iPt]->GetEntry(iEntry);
                if(chiSquareBC > maxChi2 || chiSquareBC < minChi2 || TMath::Abs(nSigmaBC-nSigmaBinCounting[iBinCount]) > 0.001 || TMath::Abs(confSigmaBC-minSigmaConf) > 0.001 || TMath::Abs(confMeanBC-minMeanConf) > 0.001)
                    continue;
                hBinCountVsTrial[iBinCount][iPt]->SetBinContent(iBin, rawYieldBC);
                hBinCountVsTrial[iBinCount][iPt]->SetBinError(iBin, rawYielBCUnc);
                iBin++;
                Counter++;
            }
        }

        double mean = hRawYield[iPt]->GetMean();
        double RMS = hRawYield[iPt]->GetRMS();
        double disp = (max(hRawYieldVsTrial[iPt]) - min(hRawYieldVsTrial[iPt])) / TMath::Sqrt(12);
        double RMSunc = RMS / hRawYieldRef->GetBinContent(iPt+1) * 100;
        double Flatunc = disp / hRawYieldRef->GetBinContent(iPt+1) * 100;

        stats[iPt] = new TPaveText(0.15, 0.65, 0.44, 0.85, "NDC");
        stats[iPt]->SetTextSize(0.05);
        stats[iPt]->SetFillColor(0);
        stats[iPt]->SetFillStyle(0);
        stats[iPt]->SetBorderSize(0);
        stats[iPt]->SetTextFont(42);
        stats[iPt]->SetTextColor(kBlue + 1);
        stats[iPt]->AddText(Form("mean = %0.1f", mean));
        stats[iPt]->AddText(Form("RMS = %0.1f (%0.1f%%)", RMS, RMSunc));
        stats[iPt]->AddText(Form("#frac{max-min}{#sqrt{12}} = %.1f (%0.1f%%)", disp, Flatunc));

        for (unsigned int iBinCount = 0; iBinCount < nBinCounting; iBinCount++) {

            mean = hBinCount[iBinCount][iPt]->GetMean();
            RMS = hBinCount[iBinCount][iPt]->GetRMS();
            disp = (max(hBinCountVsTrial[iBinCount][iPt]) - min(hBinCountVsTrial[iBinCount][iPt])) / TMath::Sqrt(12);
            RMSunc = RMS / hRawYieldRef->GetBinContent(iPt+1) * 100;
            Flatunc = disp / hRawYieldRef->GetBinContent(iPt+1) * 100;

            statsbc[iBinCount][iPt] = new TPaveText(0.6, 0.65 - 0.3 * iBinCount, 0.89, 0.85 - 0.3 * iBinCount, "NDC");
            statsbc[iBinCount][iPt]->SetTextSize(0.05);
            statsbc[iBinCount][iPt]->SetFillColor(0);
            statsbc[iBinCount][iPt]->SetFillStyle(0);
            statsbc[iBinCount][iPt]->SetBorderSize(0);
            statsbc[iBinCount][iPt]->SetTextFont(42);
            statsbc[iBinCount][iPt]->SetTextColor(binCountColors[iBinCount]);
            statsbc[iBinCount][iPt]->AddText(Form("mean = %0.1f", mean));
            statsbc[iBinCount][iPt]->AddText(Form("RMS = %0.1f (%0.1f%%)", RMS, RMSunc));
            statsbc[iBinCount][iPt]->AddText(Form("#frac{max-min}{#sqrt{12}} = %0.1f (%0.1f%%)", disp, Flatunc));
        }

        if (iPt == 0) {
            leg->AddEntry(lRawRefVsTrial[iPt], "Central value", "l");
            leg->AddEntry(hRawYieldVsTrial[iPt], "Fit method", "lpe");
            for (unsigned int iBinCount = 0; iBinCount < nBinCounting; iBinCount++) {
                leg->AddEntry(hBinCountVsTrial[iBinCount][iPt], Form("Bin counting (%0.1f#sigma)", nSigmaBinCounting[iBinCount]), "lpe");
            }
        }

        if(gPad)
            gPad->Close();
    }

    // produce plots
    TCanvas *cMultiTrial[nPtBins];
    for(int iPt = 0; iPt < nPtBins; iPt++) {
        cMultiTrial[iPt] = new TCanvas(Form("cMultiTrial_pT_%0.f-%0.f", PtLims[iPt], PtLims[iPt+1]), "", 1920, 1080);
        cMultiTrial[iPt]->Divide(2, 2);

        cMultiTrial[iPt]->cd(1);
        hBinCountVsTrial[0][iPt]->Draw();
        for (unsigned int iBinCount = 0; iBinCount < nBinCounting; iBinCount++) {
            hBinCountVsTrial[iBinCount][iPt]->Draw("same");
        }
        hRawYieldVsTrial[iPt]->Draw("same");
        if (lRawRefVsTrial[iPt]) {
            lRawRefVsTrial[iPt]->Draw("same");
        }
        leg->Draw("same");

        cMultiTrial[iPt]->cd(2);
        hBinCount[0][iPt]->GetYaxis()->SetRangeUser(0., hRawYield[iPt]->GetMaximum() * 1.5);
        hBinCount[0][iPt]->Draw();
        for (unsigned int iBinCount = 0; iBinCount < nBinCounting; iBinCount++) {
            hBinCount[iBinCount][iPt]->Draw("same");
        }
        hRawYield[iPt]->Draw("same");
        if (loadref != 1 && loadref != 2 && loadref != 4) {
            lRawRef[iPt] = new TLine(hRawYieldRef->GetBinContent(iPt+1), 0, hRawYieldRef->GetBinContent(iPt+1), hRawYield[iPt]->GetMaximum() * 1.5);
            lRawRef[iPt]->SetLineColor(kRed);
            lRawRef[iPt]->SetLineWidth(2);
            lRawRef[iPt]->Draw("same");
        }
        stats[iPt]->Draw("same");
        for (unsigned int iBinCount = 0; iBinCount < nBinCounting; iBinCount++) {
            statsbc[iBinCount][iPt]->Draw("same");
        }

        cMultiTrial[iPt]->cd(3);
        hSigmaVsTrial[iPt]->Draw();
        if (lSigmaRefVsTrial[iPt]) {
            lSigmaRefVsTrial[iPt]->Draw("same");
        }
        cMultiTrial[iPt]->cd(4);
        hChiSquareVsTrial[iPt]->Draw("p");
    }

    // save output files
    TFile outFile(outFileName.data(), "recreate");
    for (int iPt = 0; iPt < nPtBins; iPt++) {
        ntupleFit[iPt]->Write();
        ntupleBinC[iPt]->Write();
        hRawYield[iPt]->Write();
        hSigma[iPt]->Write();
        hMean[iPt]->Write();
        hChiSquare[iPt]->Write();
        hRawYieldVsTrial[iPt]->Write();
        hSigmaVsTrial[iPt]->Write();
        hMeanVsTrial[iPt]->Write();
        hChiSquareVsTrial[iPt]->Write();
        for (unsigned int iBinCount = 0; iBinCount < nBinCounting; iBinCount++) {
            hBinCount[iBinCount][iPt]->Write();
            hBinCountVsTrial[iBinCount][iPt]->Write();
        }
        cMultiTrial[iPt]->Write();
    }
    outFile.Close();
    std::cout << "\n" << outFileName.data() << " saved." << std::endl;

    outFileName = std::regex_replace(outFileName, std::regex(".root"), ".pdf");
    cMultiTrial[0]->SaveAs(Form("%s[", outFileName.data()));
    for (int iPt = 0; iPt < nPtBins; iPt++) 
        cMultiTrial[iPt]->SaveAs(outFileName.data());
    cMultiTrial[nPtBins-1]->SaveAs(Form("%s]", outFileName.data()));

    return 0;
}


//__________________________________________________________________________________________________________________
double min(TH1F *histo) {

    double min = 1.e+20;
    for (int iBin = 0; iBin < histo->GetNbinsX(); iBin++) {
        if (histo->GetBinContent(iBin + 1) > 0 && histo->GetBinContent(iBin + 1) < min)
            min = histo->GetBinContent(iBin + 1);
    }
    return min;
}


//__________________________________________________________________________________________________________________
double max(TH1F *histo) {

    double max = -1.e+20;
    for (int iBin = 0; iBin < histo->GetNbinsX(); iBin++) {
        if (histo->GetBinContent(iBin + 1) > 0 && histo->GetBinContent(iBin + 1) > max)
            max = histo->GetBinContent(iBin + 1);
    }
    return max;
}


//__________________________________________________________________________________________________________________
int LoadRefFiles(string refFileName, string refFileNameMC,
                 TH1F *&hRawYieldRef, TH1F *&hSigmaRef,
                 TH1F *&hMeanRef, TH1F *&hSigmaMC) {

    TFile *reffile = TFile::Open(refFileName.data());
    if (reffile) {
        hRawYieldRef = (TH1F *)reffile->Get("hRawYields");
        hSigmaRef = (TH1F *)reffile->Get("hRawYieldsSigma");
        hMeanRef = (TH1F *)reffile->Get("hRawYieldsMean");
        if (hRawYieldRef)
            hRawYieldRef->SetDirectory(0);
        if (hSigmaRef)
            hSigmaRef->SetDirectory(0);
        if (hMeanRef)
            hMeanRef->SetDirectory(0);
        reffile->Close();
    }

    TFile *reffileMC = TFile::Open(refFileNameMC.data());
    if (reffileMC) {
        hSigmaMC = (TH1F *)reffileMC->Get("hRawYieldsSigma");
        if (hSigmaMC)
            hSigmaMC->SetDirectory(0);
        reffileMC->Close();
    }

    if (!reffile && !reffileMC)
        return 1;
    if (!reffile)
        return 2;
    if (!reffileMC)
        return 3;
    if (!hRawYieldRef)
        return 4;
    if (!hSigmaRef)
        return 5;
    if (!hMeanRef)
        return 6;
    if (!hSigmaMC)
        return 7;

    return 0;
}


//__________________________________________________________________________________________________________________
void SetStyle() {
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadRightMargin(0.075);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetTitleOffset(1.2, "y");
    gStyle->SetTitleOffset(1.2, "x");
    gStyle->SetTitleSize(0.055, "xy");
    gStyle->SetLabelSize(0.05, "xy");    
    gStyle->SetLegendBorderSize(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    TGaxis::SetMaxDigits(3);
}
