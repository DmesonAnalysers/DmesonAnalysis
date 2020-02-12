#if !defined(__CINT__) || defined(__CLING__)

#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>

#include "AliESDtrackCuts.h"
#include "AliRDHFCutsDstoKKpi.h"

#endif

//____________________________________________________________________________________________________//
// Methods:
// 1) MakeFileForCutsDspp5TeV_TreeML --> loose cuts for 2019 analysis (non-prompt)
//____________________________________________________________________________________________________//

//__________________________________________________________________________________________
AliRDHFCutsDstoKKpi *MakeFileForCutsDspp5TeV_TreeML(bool fIsMC = false, bool fIsQA = false, bool fUseStrongPID = false, double maxPtstrongPID = -1.0)
{

    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3, 1.e10);

    const int nptbins = 2;
    float *ptbins;
    ptbins = new float[nptbins + 1];

    ptbins[0] = 1.;
    ptbins[1] = 6.;
    ptbins[2] = 50.;

    const int nvars = 20;
    float **anacutsval;
    anacutsval = new float *[nvars];
    for (int ic = 0; ic < nvars; ic++)
        anacutsval[ic] = new float[nptbins];

    /*
     Cut list                                   Rejection condition
     0 inv. mass [GeV]",                        invmassDS-massDspdg>fCutsRD
     1 pTK [GeV/c]",                            pTK<fCutsRd
     2 pTPi [GeV/c]",                           pTPi<fCutsRd
     3 d0K [cm]",                               d0K<fCutsRd
     4 d0Pi [cm]",                              d0Pi<fCutsRd
     5 dist12 [cm]",                            dist12<fCutsRd
     6 sigmavert [cm]",                         sigmavert>fCutsRd
     7 decLen [cm]",                            decLen<fCutsRD
     8 ptMax [GeV/c]",                          ptMax<fCutsRD
     9 cosThetaPoint",                          CosThetaPoint<fCutsRD
     10 Sum d0^2 (cm^2)",                       sumd0<fCutsRD
     11 dca [cm]",                              dca(i)>fCutsRD
     12 inv. mass (Mphi-MKK) [GeV]",            invmass-pdg>fCutsRD
     13 inv. mass (MKo*-MKpi) [GeV]",           invmass-pdg>fCutsRD
     14 Abs(CosineKpiPhiRFrame)^3",
     15 CosPiDsLabFrame",
     16 decLenXY [cm]"
     17 NormdecLen",
     18 NormdecLenXY [cm]",
     19 cosThetaPointXY"
    */

    for (int ipt = 0; ipt < nptbins; ipt++)
    {
        anacutsval[0][ipt]=0.25;
        anacutsval[1][ipt]=0.3;
        anacutsval[2][ipt]=0.3;
        anacutsval[3][ipt]=0.;
        anacutsval[4][ipt]=0.;
        anacutsval[5][ipt]=0.;
        anacutsval[6][ipt]=0.08;
        anacutsval[7][ipt]=0.0;
        anacutsval[8][ipt]=0.;
        anacutsval[9][ipt]=0.85;
        anacutsval[10][ipt]=0.;
        anacutsval[11][ipt]=1000.0;
        anacutsval[12][ipt]=0.015;
        anacutsval[13][ipt]=0.001;
        anacutsval[14][ipt]=0.;
        anacutsval[15][ipt]=1.;
        anacutsval[16][ipt]=0.;
        anacutsval[17][ipt]=0.;
        anacutsval[18][ipt]=0.;
        anacutsval[19][ipt]=0.85;
    }

    AliRDHFCutsDstoKKpi *analysiscuts = new AliRDHFCutsDstoKKpi();
    if (fIsQA)
        analysiscuts->SetName("DstoKKpiCuts");
    else
        analysiscuts->SetName("AnalysisCuts");
    
    analysiscuts->SetTitle("Cuts for Ds non-prompt analysis");
    analysiscuts->SetPtBins(nptbins + 1, ptbins);
    analysiscuts->SetCuts(nvars, nptbins, anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);

    analysiscuts->SetUsePID(true);
    if (fUseStrongPID)
    {
        analysiscuts->SetPidOption(AliRDHFCutsDstoKKpi::kStrong);
        analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
    }
    else
    {
       analysiscuts->SetPidOption(AliRDHFCutsDstoKKpi::kConservative);
    }

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");
    analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    if (fIsMC)
        analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
    analysiscuts->SetRemoveDaughtersFromPrim(true);
    analysiscuts->SetMinVtxContr(1);

    analysiscuts->SetMinPtCandidate(1.);
    analysiscuts->SetMaxPtCandidate(50.);

    std::cout << "This is the object I'm going to save:" << nptbins << std::endl;

    TString triggername = "kINT7";
    if (fIsMC)
        triggername = "kMB";
    TString PIDsuffix = "";
    if (fUseStrongPID)
        PIDsuffix = Form("_strongPIDpt%0.f", maxPtstrongPID);
    TString QAsuffix = "";
    if (fIsQA)
        QAsuffix = "_forQA";

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = Form("DstoKKpiCuts_pp_nonprompt_loose%s_%s%s.root", PIDsuffix.Data(), triggername.Data(), QAsuffix.Data());
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}
