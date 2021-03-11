#if !defined(__CINT__) || defined(__CLING__)

#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>

#include "AliESDtrackCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"

#endif

//____________________________________________________________________________________________________//
// Methods:
// 1) MakeFileForCutsDpluspp5TeV_TreeML --> loose cuts for 2020 analysis (femto in 1 < pT < 10 GeV/c)
// 2) MakeFileForCutsDpluspp5TeV_TreeML_LowPt --> loose cuts for 2020 analysis (femto in 0 < pT < 1 GeV/c)
// 3) MakeFileForCutsDpluspp5TeV_TreeML_AllPt --> loose cuts for 2020 analysis (femto in 0 < pT < 10 GeV/c)
// 4) MakeFileForCutsDpluspp5TeV_BkgCuts --> loose cuts for background studies
//____________________________________________________________________________________________________//

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi *MakeFileForCutsDpluspp5TeV_TreeML(bool fIsMC = false)
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
    ptbins[1] = 5.;
    ptbins[2] = 10.;

    const int nvars = 14;
    float **anacutsval;
    anacutsval = new float *[nvars];
    for (int ic = 0; ic < nvars; ic++)
        anacutsval[ic] = new float[nptbins];

    /*
     Cut list
     0          "inv. mass [GeV]",
     1			"pTK [GeV/c]",
     2			"pTPi [GeV/c]",
     3			"d0K [cm]",
     4			"d0Pi [cm]",
     5			"dist12 [cm]",
     6			"sigmavert [cm]",
     7			"decLen [cm]",
     8			"ptMax [GeV/c]",
     9			"cosThetaPoint",
     10			"Sum d0^2 (cm^2)",
     11			"dca [cm]",
     12			"norm decay length XY",
     13			"cosThetaPointXY";
     */

    for (int ipt = 0; ipt < nptbins; ipt++)
    {
        anacutsval[0][ipt]  = 0.2;   //minv
        anacutsval[1][ipt]  = 0.3;   //ptK
        anacutsval[2][ipt]  = 0.3;   //ptPi
        anacutsval[3][ipt]  = 0.0;   //d0K
        anacutsval[4][ipt]  = 0.0;   //d0Pi
        anacutsval[5][ipt]  = 0.0;   //dist12
        anacutsval[8][ipt]  = 0.0;   //pM
        anacutsval[10][ipt] = 0.0;   //sumd02
        anacutsval[11][ipt] = 1.e10; //dca
        anacutsval[12][ipt] = 0.;   //ndlXY
    }

    //pT 1-5
    anacutsval[6][0] = 0.040; //sigvert
    anacutsval[7][0] = 0.030; //declen
    anacutsval[9][0] = 0.85;  //cosp
    anacutsval[13][0] = 0.90; //cospXY

    //pT 5-50
    anacutsval[6][1] = 0.060; //sigvert
    anacutsval[7][1] = 0.030; //declen
    anacutsval[9][1] = 0.75;  //cosp
    anacutsval[13][1] = 0.80; //cospXY

    AliRDHFCutsDplustoKpipi *analysiscuts = new AliRDHFCutsDplustoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
    analysiscuts->SetPtBins(nptbins + 1, ptbins);
    analysiscuts->SetCuts(nvars, nptbins, anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetScaleNormDLxyBypOverPt(false);

    analysiscuts->SetUsePID(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");
    analysiscuts->SetTriggerMask(AliVEvent::kHighMultV0); // high multiplicity V0M
    if (fIsMC)
        analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUseImpParProdCorrCut(false);

    analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
    analysiscuts->SetRemoveDaughtersFromPrim(true);
    analysiscuts->SetMinVtxContr(1);

    analysiscuts->SetMinPtCandidate(1.);
    analysiscuts->SetMaxPtCandidate(10.);

    std::cout << "This is the object I'm going to save:" << nptbins << std::endl;

    TString triggername = "kHighMultV0";
    if (fIsMC)
        triggername = "kMB";
    TString PIDsuffix = "";

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = Form("DplustoKpipiCuts_pp_femto_loose%s_%s.root", PIDsuffix.Data(), triggername.Data());
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi *MakeFileForCutsDpluspp5TeV_TreeML_LowPt(bool fIsMC = false)
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

    ptbins[0] = 0.;
    ptbins[1] = 1.;

    const int nvars = 14;
    float **anacutsval;
    anacutsval = new float *[nvars];
    for (int ic = 0; ic < nvars; ic++)
        anacutsval[ic] = new float[nptbins];

    //pT 0-1
    anacutsval[0][0] = 0.200;   //minv
    anacutsval[1][0] = 0.300;   //ptK
    anacutsval[2][0] = 0.300;   //ptPi
    anacutsval[3][0] = 0.000;   //d0K
    anacutsval[4][0] = 0.000;   //d0Pi
    anacutsval[5][0] = 0.000;   //dist12
    anacutsval[6][0] = 0.040; //sigvert
    anacutsval[7][0] = 0.020; //declen
    anacutsval[8][0] = 0.000;   //pM
    anacutsval[9][0] = 0.850;  //cosp
    anacutsval[10][0] = 0.00;   //sumd02
    anacutsval[11][0] = 1.e10; //dca
    anacutsval[12][0] = 0.00;   //ndlXY
    anacutsval[13][0] = 0.90; //cospXY

    AliRDHFCutsDplustoKpipi *analysiscuts = new AliRDHFCutsDplustoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
    analysiscuts->SetPtBins(nptbins + 1, ptbins);
    analysiscuts->SetCuts(nvars, nptbins, anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetScaleNormDLxyBypOverPt(false);

    analysiscuts->SetUsePID(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");
    analysiscuts->SetTriggerMask(AliVEvent::kHighMultV0); // high multiplicity V0M
    if (fIsMC)
        analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUseImpParProdCorrCut(false);

    analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
    analysiscuts->SetRemoveDaughtersFromPrim(true);
    analysiscuts->SetMinVtxContr(1);

    analysiscuts->SetMinPtCandidate(0.);
    analysiscuts->SetMaxPtCandidate(1.);

    std::cout << "This is the object I'm going to save:" << nptbins << std::endl;

    TString triggername = "kHighMultV0";
    if (fIsMC)
        triggername = "kMB";
    TString PIDsuffix = "";

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = Form("DplustoKpipiCuts_pp_femto_loose_lowpt%s_%s.root", PIDsuffix.Data(), triggername.Data());
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi *MakeFileForCutsDpluspp5TeV_TreeML_AllPt(bool fIsMC = false)
{

    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3, 1.e10);

    const int nptbins = 3;
    float *ptbins;
    ptbins = new float[nptbins + 1];

    ptbins[0] = 0.;
    ptbins[1] = 1.;
    ptbins[2] = 5.;
    ptbins[3] = 10.;

    const int nvars = 14;
    float **anacutsval;
    anacutsval = new float *[nvars];
    for (int ic = 0; ic < nvars; ic++)
        anacutsval[ic] = new float[nptbins];

    for (int ipt = 0; ipt < nptbins; ipt++)
    {
        anacutsval[0][ipt]  = 0.2;   //minv
        anacutsval[1][ipt]  = 0.3;   //ptK
        anacutsval[2][ipt]  = 0.3;   //ptPi
        anacutsval[3][ipt]  = 0.0;   //d0K
        anacutsval[4][ipt]  = 0.0;   //d0Pi
        anacutsval[5][ipt]  = 0.0;   //dist12
        anacutsval[8][ipt]  = 0.0;   //pM
        anacutsval[10][ipt] = 0.0;   //sumd02
        anacutsval[11][ipt] = 1.e10; //dca
        anacutsval[12][ipt] = 0.;   //ndlXY
    }

    //pT 0-1
    anacutsval[6][0] = 0.040; //sigvert
    anacutsval[7][0] = 0.020; //declen
    anacutsval[9][0] = 0.850;  //cosp
    anacutsval[13][0] = 0.90; //cospXY

    //pT 1-5
    anacutsval[6][0] = 0.040; //sigvert
    anacutsval[7][0] = 0.030; //declen
    anacutsval[9][0] = 0.85;  //cosp
    anacutsval[13][0] = 0.90; //cospXY

    //pT 5-50
    anacutsval[6][1] = 0.060; //sigvert
    anacutsval[7][1] = 0.030; //declen
    anacutsval[9][1] = 0.75;  //cosp
    anacutsval[13][1] = 0.80; //cospXY

    AliRDHFCutsDplustoKpipi *analysiscuts = new AliRDHFCutsDplustoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
    analysiscuts->SetPtBins(nptbins + 1, ptbins);
    analysiscuts->SetCuts(nvars, nptbins, anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetScaleNormDLxyBypOverPt(false);

    analysiscuts->SetUsePID(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");
    analysiscuts->SetTriggerMask(AliVEvent::kHighMultV0); // high multiplicity V0M
    if (fIsMC)
        analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUseImpParProdCorrCut(false);

    analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
    analysiscuts->SetRemoveDaughtersFromPrim(true);
    analysiscuts->SetMinVtxContr(1);

    analysiscuts->SetMinPtCandidate(0.);
    analysiscuts->SetMaxPtCandidate(10.);

    std::cout << "This is the object I'm going to save:" << nptbins << std::endl;

    TString triggername = "kHighMultV0";
    if (fIsMC)
        triggername = "kMB";
    TString PIDsuffix = "";

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = Form("DplustoKpipiCuts_pp_femto_loose_allpt%s_%s.root", PIDsuffix.Data(), triggername.Data());
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi *MakeFileForCutsDpluspp5TeV_Bkg(bool fIsMC = false, int cutOpt = 0)
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
    ptbins[1] = 5.;
    ptbins[2] = 10.;

    const int nvars = 14;
    float **anacutsval;
    anacutsval = new float *[nvars];
    for (int ic = 0; ic < nvars; ic++)
        anacutsval[ic] = new float[nptbins];

    /*
     Cut list
     0          "inv. mass [GeV]",
     1			"pTK [GeV/c]",
     2			"pTPi [GeV/c]",
     3			"d0K [cm]",
     4			"d0Pi [cm]",
     5			"dist12 [cm]",
     6			"sigmavert [cm]",
     7			"decLen [cm]",
     8			"ptMax [GeV/c]",
     9			"cosThetaPoint",
     10			"Sum d0^2 (cm^2)",
     11			"dca [cm]",
     12			"norm decay length XY",
     13			"cosThetaPointXY";
     */

    for (int ipt = 0; ipt < nptbins; ipt++)
    {
        anacutsval[0][ipt]  = 0.2;   //minv
        anacutsval[1][ipt]  = 0.3;   //ptK
        anacutsval[2][ipt]  = 0.3;   //ptPi
        anacutsval[3][ipt]  = 0.0;   //d0K
        anacutsval[4][ipt]  = 0.0;   //d0Pi
        anacutsval[5][ipt]  = 0.0;   //dist12
        anacutsval[6][ipt]  = 0.04;  //sigmavert
        anacutsval[7][ipt]  = 0.04;  //decLen
        anacutsval[8][ipt]  = 0.0;   //pM
        anacutsval[9][ipt]  = 0.95;  //cosp
        anacutsval[10][ipt] = 0.0;   //sumd02
        anacutsval[11][ipt] = 1.e10; //dca
        anacutsval[12][ipt] = 0.;    //ndlXY
    }

    if(cutOpt == 0)
    {
        //pT 1-5
        anacutsval[12][0] = 4.;   //ndlXY
        anacutsval[13][0] = 0.95; //cospXY
        //pT 5-10
        anacutsval[12][1] = 4.;   //ndlXY
        anacutsval[13][1] = 0.95; //cospXY
    }
    else if(cutOpt == 1)
    {
        //pT 1-5
        anacutsval[12][0] = 6.;   //ndlXY
        anacutsval[13][0] = 0.97; //cospXY
        //pT 5-10
        anacutsval[12][1] = 6.;   //ndlXY
        anacutsval[13][1] = 0.97; //cospXY
    }
    else if(cutOpt == 2)
    {
        //pT 1-5
        anacutsval[12][0] = 8.;   //ndlXY
        anacutsval[13][0] = 0.99; //cospXY
        //pT 5-10
        anacutsval[12][1] = 8.;   //ndlXY
        anacutsval[13][1] = 0.99; //cospXY
    }
    else if(cutOpt == 3)
    {
        //pT 1-5
        anacutsval[12][0] = 10.;   //ndlXY
        anacutsval[13][0] = 0.995; //cospXY
        //pT 5-10
        anacutsval[12][1] = 10.;   //ndlXY
        anacutsval[13][1] = 0.995; //cospXY
    }
    
    AliRDHFCutsDplustoKpipi *analysiscuts = new AliRDHFCutsDplustoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
    analysiscuts->SetPtBins(nptbins + 1, ptbins);
    analysiscuts->SetCuts(nvars, nptbins, anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetScaleNormDLxyBypOverPt(false);

    analysiscuts->SetUsePID(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");
    analysiscuts->SetTriggerMask(AliVEvent::kHighMultV0); // high multiplicity V0M
    if (fIsMC)
        analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUseImpParProdCorrCut(false);

    analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
    analysiscuts->SetRemoveDaughtersFromPrim(true);
    analysiscuts->SetMinVtxContr(1);

    analysiscuts->SetMinPtCandidate(1.);
    analysiscuts->SetMaxPtCandidate(50.);

    std::cout << "This is the object I'm going to save:" << nptbins << std::endl;

    TString triggername = "kHighMultV0";
    if (fIsMC)
        triggername = "kMB";
    TString PIDsuffix = "";

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = Form("DplustoKpipiCuts_pp_femto_bkgStudy%d%s_%s.root", cutOpt, PIDsuffix.Data(), triggername.Data());
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}
