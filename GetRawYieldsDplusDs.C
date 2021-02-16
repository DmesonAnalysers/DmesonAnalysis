//___________________________________________________________________________________//
// Macro for fitting D+ and Ds+ invariant-mass spectra                               //
// Main Function: GetRawYieldsDplusDs                                                //
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

enum {k010, k3050, k6080, kpp5TeVPrompt, kpp5TeVFD, kpp13TeVPrompt, kpp13TeVFD};
enum {kDplus, kDs, kLc};

//__________________________________________________________________________________________________________________
int GetRawYieldsDplusDs(int cent = k010, bool isMC = false, TString infilename = "InvMassSpectraDplus_010_PbPb2015cuts.root", TString cfgfilename = "Dplus/config_Dplus_Fit.yml", TString outFileName = "RawYieldsDplus_010_PbPb2015cuts.root");
double SingleGaus(double *m, double *pars);
double DoublePeakSingleGaus(double *x, double *pars);
double DoubleGaus(double *m, double *pars);
double DoublePeakDoubleGaus(double *m, double *pars);
void SetHistoStyle(TH1 *histo, int color=kBlack, double markersize=1.);
void SetStyle();
void DivideCanvas(TCanvas* c, int nPtBins);

//__________________________________________________________________________________________________________________
int GetRawYieldsDplusDs(int cent, bool isMC, TString infilename, TString cfgfilename, TString outFileName) {
    SetStyle();

    //load config
    TString centname = "";
    if(cent==k010) centname = "Cent010";
    else if(cent==k3050) centname = "Cent3050";
    else if(cent==k6080) centname = "Cent6080";
    else if(cent==kpp5TeVPrompt) centname = "pp5TeVPrompt";
    else if(cent==kpp5TeVFD) centname = "pp5TeVFD";
    else if(cent==kpp13TeVPrompt) centname = "pp13TeVPrompt";
    else if(cent==kpp13TeVFD) centname = "pp13TeVFD";
    else {
        cerr << "ERROR: centrality "<< cent << " is not supported! Exit"<<endl;
        return -1;
    } 
    YAML::Node config = YAML::LoadFile(cfgfilename.Data());

    string ParticleName = config[centname.Data()]["Particle"].as<string>();
    int particle;
    if(ParticleName=="Dplus"){
        particle=kDplus;
    } else if (ParticleName=="Ds"){
        particle=kDs;
    } else if (ParticleName=="Lc"){
        particle=kLc;
    } else{
        cerr << "ERROR: only Dplus, Ds and Lc are supported! Exit";
        return -1;
    }
    
    bool fixSigma = static_cast<bool>(config[centname.Data()]["FixSigma"].as<int>());
    string infilenameSigma = config[centname.Data()]["SigmaFile"].as<string>();
    bool isSigmaMultFromUnc = false;
    double sigmaMult = 1.;
    string sigmaMultFromUnc = "";
    try {
        sigmaMult = config[centname.Data()]["SigmaMultFactor"].as<double>();
    } catch (const YAML::BadConversion& e)  {
        sigmaMultFromUnc = config[centname.Data()]["SigmaMultFactor"].as<string>();
        isSigmaMultFromUnc = true;
    }
    bool fixMean = static_cast<bool>(config[centname.Data()]["FixMean"].as<int>());
    string infilenameMean = config[centname.Data()]["MeanFile"].as<string>();
    bool UseLikelihood = config[centname.Data()]["UseLikelihood"].as<int>();
    vector<double> PtMin = config[centname.Data()]["PtMin"].as<vector<double>>();
    vector<double> PtMax = config[centname.Data()]["PtMax"].as<vector<double>>();
    vector<double> MassMin = config[centname.Data()]["MassMin"].as<vector<double>>();
    vector<double> MassMax = config[centname.Data()]["MassMax"].as<vector<double>>();
    vector<int> Rebin = config[centname.Data()]["Rebin"].as<vector<int>>();
    vector<int> InclSecPeak = config[centname.Data()]["InclSecPeak"].as<vector<int>>();
    vector<double> SigmaSecPeak = config[centname.Data()]["SigmaSecPeak"].as<vector<double>>();
    string infilenameSigmaSecPeak = config[centname.Data()]["SigmaFileSecPeak"].as<string>();
    double sigmaMultSecPeak = config[centname.Data()]["SigmaMultFactorSecPeak"].as<double>();
    bool fixSigmaToFirstPeak = static_cast<bool>(config[centname.Data()]["FixSigmaToFirstPeak"].as<int>());
    vector<string> bkgfunc = config[centname.Data()]["BkgFunc"].as<vector<string>>();
    vector<string> sgnfunc = config[centname.Data()]["SgnFunc"].as<vector<string>>();
    const unsigned int nPtBins = PtMin.size();
    int BkgFunc[nPtBins], SgnFunc[nPtBins], degPol[nPtBins];
    double PtLims[nPtBins+1];

    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
        PtLims[iPt] = PtMin[iPt];
        PtLims[iPt+1] = PtMax[iPt];

        degPol[iPt] = -1;
        if(bkgfunc[iPt] == "kExpo")
            BkgFunc[iPt] = AliHFInvMassFitter::kExpo;
        else if(bkgfunc[iPt] == "kLin")
            BkgFunc[iPt] = AliHFInvMassFitter::kLin;
        else if(bkgfunc[iPt] == "kPol2")
            BkgFunc[iPt] = AliHFInvMassFitter::kPol2;
        else if(bkgfunc[iPt] == "kPol3")
        {
            BkgFunc[iPt] = 6;
            degPol[iPt] = 3;
            if (PtMin.size() > 1 && InclSecPeak[iPt] == 1) {
                cerr << "Pol3 and Pol4 fits work only with one bin if you have the secondary peak! Exit!" << endl;
                return -1;
            }
        }
        else if(bkgfunc[iPt] == "kPol4")
        {
            BkgFunc[iPt] = 6;
            degPol[iPt] = 4;
            if (PtMin.size() > 1 && InclSecPeak[iPt] == 1) {
                cerr << "Pol3 and Pol4 fits work only with one bin if you have the secondary peak! Exit!" << endl;
                return -1;
            }
        }
        else
        {
            cerr << "ERROR: only kExpo, kLin, kPol2, kPol3, and kPol4 background functions supported! Exit" << endl;
            return -1;
        }
        

        if(sgnfunc[iPt] == "kGaus")
            SgnFunc[iPt] = AliHFInvMassFitter::kGaus;
        else if(sgnfunc[iPt] == "k2Gaus")
            SgnFunc[iPt] = AliHFInvMassFitter::k2Gaus;
        else
        {
            cerr << "ERROR: only kGaus and k2Gaus signal functions supported! Exit" << endl;
            return -1;
        }
    }

    TString massaxistit = "";
    if(particle==kDplus) massaxistit = "#it{M}(K#pi#pi) (GeV/#it{c}^{2})";
    else if(particle==kDs) massaxistit = "#it{M}(KK#pi) (GeV/#it{c}^{2})";
    else if(particle==kLc) massaxistit = "#it{M}(pK^{0}_{s}) (GeV/#it{c}^{2})";
    
    //load inv-mass histos
    auto infile = TFile::Open(infilename.Data());
    if(!infile || !infile->IsOpen()) return -1;
    TH1F* hMass[nPtBins];
    TH1F* hEv = NULL;

    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {
        if(!isMC)
            hMass[iPt] = static_cast<TH1F*>(infile->Get(Form("hMass_%0.f_%0.f",PtMin[iPt]*10,PtMax[iPt]*10)));
        else {
            hMass[iPt] = static_cast<TH1F*>(infile->Get(Form("hPromptMass_%0.f_%0.f",PtMin[iPt]*10,PtMax[iPt]*10)));
            hMass[iPt]->Add(static_cast<TH1F*>(infile->Get(Form("hFDMass_%0.f_%0.f",PtMin[iPt]*10,PtMax[iPt]*10))));
            if(InclSecPeak[iPt]) {
                hMass[iPt]->Add(static_cast<TH1F*>(infile->Get(Form("hPromptSecPeakMass_%0.f_%0.f",PtMin[iPt]*10,PtMax[iPt]*10))));
                hMass[iPt]->Add(static_cast<TH1F*>(infile->Get(Form("hFDSecPeakMass_%0.f_%0.f",PtMin[iPt]*10,PtMax[iPt]*10))));
            }
        }
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
    auto hRawYieldsSignificanceSecondPeak = new TH1D("hRawYieldsSignificanceSecondPeak",";#it{p}_{T} (GeV/#it{c});signficance second peak (3#sigma)",nPtBins,PtLims);
    auto hRawYieldsSigmaRatioSecondFirstPeak = new TH1D("hRawYieldsSigmaRatioSecondFirstPeak",";#it{p}_{T} (GeV/#it{c});width second peak / width first peak",nPtBins,PtLims);
    auto hRawYieldsSoverBSecondPeak = new TH1D("hRawYieldsSoverBSecondPeak",";#it{p}_{T} (GeV/#it{c});S/B second peak (3#sigma)",nPtBins,PtLims);
    auto hRawYieldsSignalSecondPeak = new TH1D("hRawYieldsSignalSecondPeak",";#it{p}_{T} (GeV/#it{c});Signal second peak (3#sigma)",nPtBins,PtLims);
    auto hRawYieldsBkgSecondPeak = new TH1D("hRawYieldsBkgSecondPeak",";#it{p}_{T} (GeV/#it{c});Background second peak (3#sigma)",nPtBins,PtLims);
    auto hRawYieldsTrue = new TH1D("hRawYieldsTrue",";#it{p}_{T} (GeV/#it{c});true signal",nPtBins,PtLims);
    auto hRawYieldsSecondPeakTrue = new TH1D("hRawYieldsSecondPeakTrue",";#it{p}_{T} (GeV/#it{c});true signal second peak",nPtBins,PtLims);
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
    SetHistoStyle(hRawYieldsSignificanceSecondPeak,kRed+1);
    SetHistoStyle(hRawYieldsSigmaRatioSecondFirstPeak,kRed+1);
    SetHistoStyle(hRawYieldsSoverBSecondPeak,kRed+1);
    SetHistoStyle(hRawYieldsSignalSecondPeak,kRed+1);
    SetHistoStyle(hRawYieldsBkgSecondPeak,kRed+1);
    SetHistoStyle(hRawYieldsTrue);
    SetHistoStyle(hRawYieldsSecondPeakTrue,kRed+1);
    SetHistoStyle(hRelDiffRawYieldsFitTrue);
    SetHistoStyle(hRelDiffRawYieldsSecondPeakFitTrue,kRed+1);

    // additional S, B, S/B, and significance histos for different Nsigma values (filled only in case of data)
    const int nMassWindows = 6;
    double nSigma4SandB[nMassWindows] = {1.5, 1.75, 2., 2.25, 2.5, 2.75};
    TH1D *hRawYieldsSignalDiffSigma[nMassWindows], *hRawYieldsBkgDiffSigma[nMassWindows];
    TH1D *hRawYieldsSoverBDiffSigma[nMassWindows], *hRawYieldsSignifDiffSigma[nMassWindows];
    for(int iS = 0; iS < nMassWindows; iS++)
    {
        hRawYieldsSignalDiffSigma[iS] = new TH1D(Form("hRawYieldsSignal_%0.2fsigma", nSigma4SandB[iS]),
                                                 Form(";#it{p}_{T} (GeV/#it{c});Signal (%0.2f#sigma)", nSigma4SandB[iS]), nPtBins, PtLims);
        hRawYieldsBkgDiffSigma[iS] = new TH1D(Form("hRawYieldsBkg_%0.2fsigma", nSigma4SandB[iS]),
                                              Form(";#it{p}_{T} (GeV/#it{c});Background (%0.2f#sigma)", nSigma4SandB[iS]), nPtBins, PtLims);
        hRawYieldsSoverBDiffSigma[iS] = new TH1D(Form("hRawYieldsSoverB_%0.2fsigma", nSigma4SandB[iS]),
                                                 Form(";#it{p}_{T} (GeV/#it{c});S/B (%0.2f#sigma)", nSigma4SandB[iS]), nPtBins, PtLims);
        hRawYieldsSignifDiffSigma[iS] = new TH1D(Form("hRawYieldsSignif_%0.2fsigma", nSigma4SandB[iS]),
                                                 Form(";#it{p}_{T} (GeV/#it{c});significance (%0.2f#sigma)", nSigma4SandB[iS]), nPtBins, PtLims);
        SetHistoStyle(hRawYieldsSignalDiffSigma[iS]);
        SetHistoStyle(hRawYieldsBkgDiffSigma[iS]);
        SetHistoStyle(hRawYieldsSoverBDiffSigma[iS]);
        SetHistoStyle(hRawYieldsSignifDiffSigma[iS]);
    }

    TH1D *hSigmaToFix = NULL;
    if(fixSigma) {
        auto infileSigma = TFile::Open(infilenameSigma.data());
        if(!infileSigma)
            return -2;
        hSigmaToFix = static_cast<TH1D*>(infileSigma->Get("hRawYieldsSigma"));
        hSigmaToFix->SetDirectory(0);
        if(static_cast<unsigned int>(hSigmaToFix->GetNbinsX()) != nPtBins)
            cout << "WARNING: Different number of bins for this analysis and histo for fix sigma" << endl;
        infileSigma->Close();
    }

    TH1D *hMeanToFix = NULL;
    if(fixMean) {
        auto infileMean = TFile::Open(infilenameMean.data());
        if(!infileMean)
            return -3;
        hMeanToFix = static_cast<TH1D*>(infileMean->Get("hRawYieldsMean"));
        hMeanToFix->SetDirectory(0);
        if(static_cast<unsigned int>(hMeanToFix->GetNbinsX())!=nPtBins)
            cout << "WARNING: Different number of bins for this analysis and histo for fix mean" << endl;
        infileMean->Close();
    }

    TH1D *hSigmaFirstPeakMC = NULL;
    TH1D *hSigmaToFixSecPeak = NULL;
    auto infileSigmaSecPeak = TFile::Open(infilenameSigmaSecPeak.data());
    if(!infileSigmaSecPeak && fixSigmaToFirstPeak)
        return -2;
    if(infileSigmaSecPeak) {
        hSigmaFirstPeakMC = static_cast<TH1D*>(infileSigmaSecPeak->Get("hRawYieldsSigma"));
        hSigmaToFixSecPeak = static_cast<TH1D*>(infileSigmaSecPeak->Get("hRawYieldsSigmaSecondPeak"));
        hSigmaFirstPeakMC->SetDirectory(0);
        hSigmaToFixSecPeak->SetDirectory(0);
        if(static_cast<unsigned int>(hSigmaFirstPeakMC->GetNbinsX()) != nPtBins ||
           static_cast<unsigned int>(hSigmaToFixSecPeak->GetNbinsX()) != nPtBins)
            cout << "WARNING: Different number of bins for this analysis and histos for fix sigma" << endl;
        infileSigmaSecPeak->Close();
    }

    //fit histos
    double massDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    double massDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
    double massLc = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
    double massForFit = -1.;
    
    if(particle==kDplus) massForFit = massDplus;
    else if (particle==kDs) massForFit = massDs;
    else if (particle==kLc) massForFit = massLc;

    TH1F* hMassForFit[nPtBins];

    int canvSize[2] = {1920, 1080};
    if(nPtBins == 1)
    {
        canvSize[0] = 500;
        canvSize[1] = 500;
    }

    int nMaxCanvases = 20; // do not put more than 20 bins per canvas to make them visible
    const int nCanvases = ceil((float)nPtBins / nMaxCanvases);
    TCanvas *cMass[nCanvases], *cResiduals[nCanvases];
    for(int iCanv=0; iCanv<nCanvases; iCanv++)
    {
        int nPads = (nCanvases == 1) ? nPtBins : nMaxCanvases;
        cMass[iCanv] = new TCanvas(Form("cMass%d", iCanv), Form("cMass%d", iCanv),canvSize[0], canvSize[1]);
        DivideCanvas(cMass[iCanv], nPads);
        cResiduals[iCanv] = new TCanvas(Form("cResiduals%d", iCanv), Form("cResiduals%d", iCanv), canvSize[0], canvSize[1]);
        DivideCanvas(cResiduals[iCanv], nPads);
    }

    for(unsigned int iPt=0; iPt<nPtBins; iPt++) {

        int iCanv = floor((float)iPt / nMaxCanvases);

        hMassForFit[iPt]=reinterpret_cast<TH1F*>(AliVertexingHFUtils::RebinHisto(hMass[iPt],Rebin[iPt]));
        TString pttitle = Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c}",PtMin[iPt],PtMax[iPt]);
        hMassForFit[iPt]->SetTitle(Form("%s;%s;Counts per %0.f MeV/#it{c}^{2}",pttitle.Data(),massaxistit.Data(),hMassForFit[iPt]->GetBinWidth(1)*1000));
        hMassForFit[iPt]->SetName(Form("MassForFit%d",iPt));
        double markerSize = 1.;
        if(nPtBins>15)
            markerSize = 0.5;
        SetHistoStyle(hMassForFit[iPt], kBlack, markerSize);

        if(isMC) { //MC
            int parRawYield = 0, parMean = 1., parSigma1 = 2; //always the same
            int parSigma2 = -1, parFrac2Gaus = -1, parRawYieldSecPeak = -1, parMeanSecPeak = -1, parSigmaSecPeak = -1;
            TF1* massFunc = NULL;
            if(SgnFunc[iPt]==AliHFInvMassFitter::kGaus) {
                if(!(InclSecPeak[iPt] && particle==kDs)) {
                    massFunc = new TF1(Form("massFunc%d",iPt),SingleGaus,MassMin[iPt],MassMax[iPt],3);
                    massFunc->SetParameters(hMassForFit[iPt]->Integral()*hMassForFit[iPt]->GetBinWidth(1),massForFit,0.010);
                }
                else {
                    massFunc = new TF1(Form("massFunc%d",iPt),DoublePeakSingleGaus,MassMin[iPt],MassMax[iPt],6);
                    massFunc->SetParameters(hMassForFit[iPt]->Integral()*hMassForFit[iPt]->GetBinWidth(1),massForFit,0.010,hMassForFit[iPt]->Integral()*hMassForFit[iPt]->GetBinWidth(1),massDplus,0.010);
                    parRawYieldSecPeak = 3;
                    parMeanSecPeak = 4;
                    parSigmaSecPeak = 5;
                }
            }
            else if(SgnFunc[iPt]==AliHFInvMassFitter::k2Gaus){
                parSigma2 = 3;
                parFrac2Gaus = 4;
                if(!(InclSecPeak[iPt] && particle==kDs)) {
                    massFunc = new TF1(Form("massFunc%d",iPt),DoubleGaus,MassMin[iPt],MassMax[iPt],5);
                    massFunc->SetParameters(hMassForFit[iPt]->Integral()*hMassForFit[iPt]->GetBinWidth(1),massForFit,0.010,0.030,0.9);
                }
                else {
                    massFunc = new TF1(Form("massFunc%d",iPt),DoublePeakDoubleGaus,MassMin[iPt],MassMax[iPt],8);
                    massFunc->SetParameters(hMassForFit[iPt]->Integral()*hMassForFit[iPt]->GetBinWidth(1),massForFit,0.010,0.030,0.9,hMassForFit[iPt]->Integral()*hMassForFit[iPt]->GetBinWidth(1),massDplus,0.010);
                    parRawYieldSecPeak = 5;
                    parMeanSecPeak = 6;
                    parSigmaSecPeak = 7;
                }
            }

            if(nPtBins>1)
                cMass[iCanv]->cd(iPt-nMaxCanvases*iCanv+1);
            else
                cMass[iCanv]->cd();
            hMassForFit[iPt]->Fit(massFunc,"E"); //fit with chi2

            double rawyield = massFunc->GetParameter(parRawYield);
            double rawyielderr = massFunc->GetParError(parRawYield);
            double sigma = massFunc->GetParameter(parSigma1);
            double sigmaerr = massFunc->GetParError(parSigma1);
            double mean = massFunc->GetParameter(parMean);
            double meanerr = massFunc->GetParError(parMean);
            double redchi2 = massFunc->GetChisquare() / massFunc->GetNDF();

            hRawYields->SetBinContent(iPt+1,rawyield);
            hRawYields->SetBinError(iPt+1,rawyielderr);
            hRawYieldsSigma->SetBinContent(iPt+1,sigma);
            hRawYieldsSigma->SetBinError(iPt+1,sigmaerr);
            hRawYieldsMean->SetBinContent(iPt+1,mean);
            hRawYieldsMean->SetBinError(iPt+1,meanerr);
            hRawYieldsChiSquare->SetBinContent(iPt+1,redchi2);
            hRawYieldsChiSquare->SetBinError(iPt+1,0.);

            hRawYieldsTrue->SetBinContent(iPt+1,hMassForFit[iPt]->Integral());
            hRawYieldsTrue->SetBinError(iPt+1,TMath::Sqrt(hMassForFit[iPt]->Integral()));
            hRelDiffRawYieldsFitTrue->SetBinContent(iPt+1,rawyield-hMassForFit[iPt]->Integral());
            hRelDiffRawYieldsFitTrue->SetBinError(iPt+1,TMath::Sqrt(rawyielderr*rawyielderr+hMassForFit[iPt]->Integral()));

            if(InclSecPeak[iPt] && particle==kDs) {
                double rawyieldSecPeak = massFunc->GetParameter(parRawYieldSecPeak);
                double rawyieldSecPeakerr = massFunc->GetParError(parRawYieldSecPeak);
                double sigmasecondpeak = massFunc->GetParameter(parSigmaSecPeak);
                double sigmasecondpeakerr = massFunc->GetParError(parSigmaSecPeak);
                double meansecondpeak = massFunc->GetParameter(parMeanSecPeak);
                double meansecondpeakerr = massFunc->GetParError(parMeanSecPeak);
                hRawYieldsSecondPeak->SetBinContent(iPt+1,rawyieldSecPeak);
                hRawYieldsSecondPeak->SetBinError(iPt+1,rawyieldSecPeakerr);
                hRawYieldsMeanSecondPeak->SetBinContent(iPt+1,meansecondpeak);
                hRawYieldsMeanSecondPeak->SetBinError(iPt+1,meansecondpeakerr);
                hRawYieldsSigmaSecondPeak->SetBinContent(iPt+1,sigmasecondpeak);
                hRawYieldsSigmaSecondPeak->SetBinError(iPt+1,sigmasecondpeakerr);
                hRawYieldsSigmaRatioSecondFirstPeak->SetBinContent(iPt+1,sigmasecondpeak/sigma);
                hRawYieldsSigmaRatioSecondFirstPeak->SetBinError(iPt+1,TMath::Sqrt(sigmaerr*sigmaerr/(sigma*sigma)+sigmasecondpeakerr*sigmasecondpeakerr/(sigmasecondpeak*sigmasecondpeak))*sigmasecondpeak/sigma); //neglected correlation between parameters

                hRawYieldsSecondPeakTrue->SetBinContent(iPt+1,rawyield);
                hRelDiffRawYieldsSecondPeakFitTrue->SetBinContent(iPt+1,rawyield);
            }
            if(SgnFunc[iPt]==AliHFInvMassFitter::k2Gaus) {
                double sigma2 = massFunc->GetParameter(parSigma2);
                double sigma2err = massFunc->GetParError(parSigma2);
                double frac2gaus = massFunc->GetParameter(parFrac2Gaus);
                double frac2gauserr = massFunc->GetParError(parFrac2Gaus);
                hRawYieldsSigma2->SetBinContent(iPt+1,sigma2);
                hRawYieldsSigma2->SetBinError(iPt+1,sigma2err);
                hRawYieldsFracGaus2->SetBinContent(iPt+1,frac2gaus);
                hRawYieldsFracGaus2->SetBinError(iPt+1,frac2gauserr);
            }
        }
        else { //data
            auto massFitter = new AliHFInvMassFitter(hMassForFit[iPt] ,MassMin[iPt], MassMax[iPt], BkgFunc[iPt], SgnFunc[iPt]);
            if(degPol[iPt] > 0)
                massFitter->SetPolDegreeForBackgroundFit(degPol[iPt]);
            if(UseLikelihood)
                massFitter->SetUseLikelihoodFit();
            if(fixMean)
                massFitter->SetFixGaussianMean(hMeanToFix->GetBinContent(iPt+1));
            else
                massFitter->SetInitialGaussianMean(massForFit);
            if(fixSigma) {
                if(!isSigmaMultFromUnc)
                    massFitter->SetFixGaussianSigma(hSigmaToFix->GetBinContent(iPt+1)*sigmaMult);
                else {
                    if(sigmaMultFromUnc=="MinusUnc")
                        massFitter->SetFixGaussianSigma(hSigmaToFix->GetBinContent(iPt+1)-hSigmaToFix->GetBinError(iPt+1));
                    else if(sigmaMultFromUnc=="PlusUnc")
                        massFitter->SetFixGaussianSigma(hSigmaToFix->GetBinContent(iPt+1)+hSigmaToFix->GetBinError(iPt+1));
                    else
                        cout << "WARNING: impossible to fix sigma! Wrong mult factor set in config file!" << endl;
                }
            }
            else {
                if(hSigmaToFix)
                    massFitter->SetInitialGaussianSigma(hSigmaToFix->GetBinContent(iPt+1)*sigmaMult);
                else
                    massFitter->SetInitialGaussianSigma(0.008);
            }

            if(InclSecPeak[iPt] && particle==kDs) {
                if (hSigmaToFixSecPeak) {
                    massFitter->IncludeSecondGausPeak(massDplus, false, hSigmaToFixSecPeak->GetBinContent(iPt+1) * sigmaMultSecPeak, true);
                    if (fixSigmaToFirstPeak) {
                        // fix D+ peak to sigmaMC(D+)/sigmaMC(Ds+)*sigmaData(Ds+)
                        massFitter->MassFitter(false);
                        double sigmaFirstPeak = massFitter->GetSigma();
                        double sigmaRatioMC = hSigmaToFixSecPeak->GetBinContent(iPt+1) / hSigmaFirstPeakMC->GetBinContent(iPt+1);
                        massFitter->IncludeSecondGausPeak(massDplus, false, sigmaRatioMC * sigmaFirstPeak, true);
                    }
                } else {
                    massFitter->IncludeSecondGausPeak(massDplus, false, SigmaSecPeak[iPt], true);
                }
            }

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

            for(int iS = 0; iS < nMassWindows; iS++)
            {
                massFitter->Significance(nSigma4SandB[iS],signif,signiferr);
                massFitter->Signal(nSigma4SandB[iS],sgn,sgnerr);
                massFitter->Background(nSigma4SandB[iS],bkg,bkgerr);

                hRawYieldsSignalDiffSigma[iS]->SetBinContent(iPt+1,sgn);
                hRawYieldsSignalDiffSigma[iS]->SetBinError(iPt+1,sgnerr);
                hRawYieldsBkgDiffSigma[iS]->SetBinContent(iPt+1,bkg);
                hRawYieldsBkgDiffSigma[iS]->SetBinError(iPt+1,bkgerr);
                hRawYieldsSoverBDiffSigma[iS]->SetBinContent(iPt+1,sgn/bkg);
                hRawYieldsSoverBDiffSigma[iS]->SetBinError(iPt+1,sgn/bkg*TMath::Sqrt(sgnerr/sgn*sgnerr/sgn+bkgerr/bkg*bkgerr/bkg));
                hRawYieldsSignifDiffSigma[iS]->SetBinContent(iPt+1,signif);
                hRawYieldsSignifDiffSigma[iS]->SetBinError(iPt+1,signiferr);
            }

            TF1* fTotFunc = massFitter->GetMassFunc();
            TF1* fBkgFunc = massFitter->GetBackgroundRecalcFunc();

            double parFrac2Gaus = -1, parsecondsigma = -1;
            if(SgnFunc[iPt]==AliHFInvMassFitter::k2Gaus) {
                if(!(InclSecPeak[iPt] && particle==kDs)) {
                    parFrac2Gaus = fTotFunc->GetNpar()-2;
                    parsecondsigma = fTotFunc->GetNpar()-1;
                }
                else {
                    parFrac2Gaus = fTotFunc->GetNpar()-5;
                    parsecondsigma = fTotFunc->GetNpar()-4;
                }

                double sigma2 = fTotFunc->GetParameter(parsecondsigma);
                double sigma2err = fTotFunc->GetParError(parsecondsigma);
                double frac2gaus = fTotFunc->GetParameter(parFrac2Gaus);
                double frac2gauserr = fTotFunc->GetParError(parFrac2Gaus);
                hRawYieldsSigma2->SetBinContent(iPt+1,sigma2);
                hRawYieldsSigma2->SetBinError(iPt+1,sigma2err);
                hRawYieldsFracGaus2->SetBinContent(iPt+1,frac2gaus);
                hRawYieldsFracGaus2->SetBinError(iPt+1,frac2gauserr);
            }

            if(InclSecPeak[iPt] && particle==kDs) {
                int paryieldSecPeak = fTotFunc->GetNpar()-3;
                int parMeansecondpeak = fTotFunc->GetNpar()-2;
                int parSigmasecondpeak = fTotFunc->GetNpar()-1;

                double rawyieldSecPeak = fTotFunc->GetParameter(paryieldSecPeak)/hMassForFit[iPt]->GetBinWidth(1);
                double rawyieldSecPeakerr = fTotFunc->GetParError(paryieldSecPeak)/hMassForFit[iPt]->GetBinWidth(1);
                double meansecondpeak = fTotFunc->GetParameter(parMeansecondpeak);
                double meansecondpeakerr = fTotFunc->GetParError(parMeansecondpeak);
                double sigmasecondpeak = fTotFunc->GetParameter(parSigmasecondpeak);
                double sigmasecondpeakerr = fTotFunc->GetParError(parSigmasecondpeak);

                double bkgSecPeak = fBkgFunc->Integral(meansecondpeak-3*sigmasecondpeak,meansecondpeak+3*sigmasecondpeak)/hMassForFit[iPt]->GetBinWidth(1);
                double bkgSecPeakerr = TMath::Sqrt(bkgSecPeak);
                double signalSecPeak = fTotFunc->Integral(meansecondpeak-3*sigmasecondpeak,meansecondpeak+3*sigmasecondpeak)/hMassForFit[iPt]->GetBinWidth(1)-bkgSecPeak;
                double signalSecPeakerr = TMath::Sqrt(signalSecPeak+bkgSecPeak);
                double signifSecPeak = -1., signifSecPeakerr = -1.;
                AliVertexingHFUtils::ComputeSignificance(signalSecPeak,signalSecPeakerr,bkgSecPeak,bkgSecPeakerr,signifSecPeak,signifSecPeakerr);

                hRawYieldsSecondPeak->SetBinContent(iPt+1,rawyieldSecPeak);
                hRawYieldsSecondPeak->SetBinError(iPt+1,rawyieldSecPeakerr);
                hRawYieldsMeanSecondPeak->SetBinContent(iPt+1,meansecondpeak);
                hRawYieldsMeanSecondPeak->SetBinError(iPt+1,meansecondpeakerr);
                hRawYieldsSigmaSecondPeak->SetBinContent(iPt+1,sigmasecondpeak);
                hRawYieldsSigmaSecondPeak->SetBinError(iPt+1,sigmasecondpeakerr);
                hRawYieldsSignificanceSecondPeak->SetBinContent(iPt+1,signifSecPeak);
                hRawYieldsSignificanceSecondPeak->SetBinError(iPt+1,signifSecPeakerr);
                hRawYieldsSigmaRatioSecondFirstPeak->SetBinContent(iPt+1,sigmasecondpeak/sigma);
                hRawYieldsSigmaRatioSecondFirstPeak->SetBinError(iPt+1,TMath::Sqrt(sigmaerr*sigmaerr/(sigma*sigma)+sigmasecondpeakerr*sigmasecondpeakerr/(sigmasecondpeak*sigmasecondpeak))*sigmasecondpeak/sigma); //neglected correlation between parameters
                hRawYieldsSoverBSecondPeak->SetBinContent(iPt+1,signalSecPeak/bkgSecPeak);
                hRawYieldsSoverBSecondPeak->SetBinError(iPt+1,signalSecPeak/bkgSecPeak*TMath::Sqrt(signalSecPeakerr/signalSecPeak*signalSecPeakerr/signalSecPeak+bkgSecPeakerr/bkgSecPeak*bkgSecPeakerr/bkgSecPeak));
                hRawYieldsSignalSecondPeak->SetBinContent(iPt+1,signalSecPeak);
                hRawYieldsSignalSecondPeak->SetBinError(iPt+1,signalSecPeakerr);
                hRawYieldsBkgSecondPeak->SetBinContent(iPt+1,bkgSecPeak);
                hRawYieldsBkgSecondPeak->SetBinError(iPt+1,bkgSecPeakerr);
            }

            if(nPtBins>1)
                cMass[iCanv]->cd(iPt-nMaxCanvases*iCanv+1);
            else
                cMass[iCanv]->cd();

            hMassForFit[iPt]->GetYaxis()->SetRangeUser(hMassForFit[iPt]->GetMinimum()*0.95,hMassForFit[iPt]->GetMaximum()*1.2);
            massFitter->DrawHere(gPad);

            if(!isMC)
            {
                //residuals
                if(nPtBins>1)
                    cResiduals[iCanv]->cd(iPt-nMaxCanvases*iCanv+1);
                else
                    cResiduals[iCanv]->cd();
                massFitter->DrawHistoMinusFit(gPad);
            }
        }
        cMass[iCanv]->Modified();
        cMass[iCanv]->Update();
        cResiduals[iCanv]->Modified();
        cResiduals[iCanv]->Update();
    }

    //save output histos
    TFile outFile(outFileName.Data(),"recreate");
    for(int iCanv=0; iCanv<nCanvases; iCanv++)
    {
        cMass[iCanv]->Write();
        if(!isMC)
            cResiduals[iCanv]->Write();
    }
    for(unsigned int iPt=0; iPt<nPtBins; iPt++)
        hMass[iPt]->Write();
    hRawYields->Write();
    hRawYieldsSigma->Write();
    hRawYieldsMean->Write();
    hRawYieldsSignificance->Write();
    hRawYieldsSoverB->Write();
    hRawYieldsSignal->Write();
    hRawYieldsBkg->Write();
    hRawYieldsChiSquare->Write();
    hRawYieldsSigma2->Write();
    hRawYieldsFracGaus2->Write();
    hRawYieldsSecondPeak->Write();
    hRawYieldsMeanSecondPeak->Write();
    hRawYieldsSigmaSecondPeak->Write();
    hRawYieldsSignificanceSecondPeak->Write();
    hRawYieldsSigmaRatioSecondFirstPeak->Write();
    hRawYieldsSoverBSecondPeak->Write();
    hRawYieldsSignalSecondPeak->Write();
    hRawYieldsBkgSecondPeak->Write();
    hRawYieldsTrue->Write();
    hRawYieldsSecondPeakTrue->Write();
    hRelDiffRawYieldsFitTrue->Write();
    hRelDiffRawYieldsSecondPeakFitTrue->Write();
    hEv->Write();
    if(!isMC)
    {
        TDirectoryFile dir("SandBDiffNsigma", "SandBDiffNsigma");
        dir.Write();
        dir.cd();
        for(int iS = 0; iS < nMassWindows; iS++)
        {
            hRawYieldsSignalDiffSigma[iS]->Write();
            hRawYieldsBkgDiffSigma[iS]->Write();
            hRawYieldsSoverBDiffSigma[iS]->Write();
            hRawYieldsSignifDiffSigma[iS]->Write();
        }
        dir.Close();
    }
    outFile.Close();

    outFileName.ReplaceAll(".root",".pdf");
    TString outFileNameRes = outFileName;
    outFileNameRes.ReplaceAll(".pdf", "_Residuals.pdf");
    for(int iCanv=0; iCanv<nCanvases; iCanv++)
    {
        if(iCanv == 0 && nCanvases > 1)
            cMass[iCanv]->SaveAs(Form("%s[", outFileName.Data()));
        cMass[iCanv]->SaveAs(outFileName.Data());
        if(iCanv == nCanvases-1 && nCanvases > 1)
            cMass[iCanv]->SaveAs(Form("%s]", outFileName.Data()));

        if(!isMC)
        {
            if(iCanv == 0 && nCanvases > 1)
                cResiduals[iCanv]->SaveAs(Form("%s[", outFileNameRes.Data()));
            cResiduals[iCanv]->SaveAs(outFileNameRes.Data());
            if(iCanv == nCanvases-1 && nCanvases > 1)
                cResiduals[iCanv]->SaveAs(Form("%s]", outFileNameRes.Data()));
        }
    }

    return 0;
}

//__________________________________________________________________________________________________________________
double SingleGaus(double *m, double *pars) {
    double norm = pars[0], mean = pars[1], sigma = pars[2];

    return norm*TMath::Gaus(m[0],mean,sigma,true);
}

//__________________________________________________________________________________________________________________
double DoubleGaus(double *m, double *pars) {
    double norm = pars[0], mean = pars[1], sigma1 = pars[2], sigma1_2 = pars[3], fg = pars[4];

    return norm*((1-fg)*TMath::Gaus(m[0],mean,sigma1,true)+fg*TMath::Gaus(m[0],mean,sigma1_2,true));
}

//__________________________________________________________________________________________________________________
double DoublePeakSingleGaus(double *m, double *pars) {
    double norm1 = pars[0], mean1 = pars[1], sigma1 = pars[2]; //Ds peak
    double norm2 = pars[3], mean2 = pars[4], sigma2 = pars[5]; //Dplus peak

    return norm1*TMath::Gaus(m[0],mean1,sigma1,true) + norm2*TMath::Gaus(m[0],mean2,sigma2,true);
}

//__________________________________________________________________________________________________________________
double DoublePeakDoubleGaus(double *m, double *pars) {
    double norm1 = pars[0], mean = pars[1], sigma1 = pars[2], sigma1_2 = pars[3], fg = pars[4]; //Ds peak
    double norm2 = pars[5], mean2 = pars[6], sigma2 = pars[7]; //Dplus peak

    return norm1*((1-fg)*TMath::Gaus(m[0],mean,sigma1,true)+fg*TMath::Gaus(m[0],mean,sigma1_2,true)) + norm2*TMath::Gaus(m[0],mean2,sigma2,true);
}

//__________________________________________________________________________________________________________________
void SetHistoStyle(TH1 *histo, int color, double markersize) {
    histo->SetStats(kFALSE);
    histo->SetMarkerSize(markersize);
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
