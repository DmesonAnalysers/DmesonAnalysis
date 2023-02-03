#if !defined(__CINT__) || defined(__CLING__)

#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>

#include "AliESDtrackCuts.h"
#include "AliRDHFCutsDStartoKpipi.h"

#endif

//____________________________________________________________________________________________________//
// Methods:
// 1) MakeFileForCutsDstarpp13TeV_TreeML_Polarization --> loose cuts for 2021 analysis
// 2) MakeFileForCutsDstarpp13TeV_TreeML_Femto --> loose cuts for 2022 analysis
// 3) MakeFileForCutsDstarpp13TeV_NoDstarCuts --> loose cuts for 2023 Ds resonance analysis
//____________________________________________________________________________________________________//

//__________________________________________________________________________________________
AliRDHFCutsDStartoKpipi *MakeFileForCutsDstarpp13TeV_TreeML_Polarization(bool fIsMC = false, TString triggername="kINT7")
{

    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3, 1.e10);

    //soft pion selections
    AliESDtrackCuts* esdTrackCutsSoftPi = new AliESDtrackCuts();
    esdTrackCutsSoftPi->SetRequireSigmaToVertex(false);
    esdTrackCutsSoftPi->SetRequireTPCRefit(false);
    esdTrackCutsSoftPi->SetRequireITSRefit(false);
    esdTrackCutsSoftPi->SetMinNClustersITS(2);
    esdTrackCutsSoftPi->SetEtaRange(-0.8, +0.8);
    esdTrackCutsSoftPi->SetPtRange(0.05, 1.e10);

    const int nptbins = 2;
    float *ptbins;
    ptbins = new float[nptbins + 1];

    ptbins[0] = 3.;
    ptbins[1] = 5.;
    ptbins[2] = 50.;

    const int nvars = 16;
    float **anacutsval;
    anacutsval = new float *[nvars];
    for (int ic = 0; ic < nvars; ic++)
        anacutsval[ic] = new float[nptbins];

    /*
     Cut list
     0          "inv. mass [GeV]",
     1			"dca [cm]",
     2			"cosThetaStar",
     3			"pTK [GeV/c]",
     4			"pTPi [GeV/c]",
     5			"d0K [cm]",
     6			"d0Pi [cm]",
     7			"d0d0 [cm^2]",
     8			"cosThetaPoint",
     9			"inv. mass half width of D* [GeV]",
     10			"half width of (M_Kpipi-M_D0) [GeV]",
     11			"PtMin of pi_s [GeV/c]",
     12			"PtMax of pi_s [GeV/c]",
     13			"theta, angle between the pi_s and decay plane of the D0 [rad]",
     14			"|cosThetaPointXY|",
     15			"NormDecayLenghtXY";
     */

    for (int ipt = 0; ipt < nptbins; ipt++)
    {
        anacutsval[1][ipt]  = 0.1;   //dca
        anacutsval[2][ipt]  = 1.1;   //cost*
        anacutsval[3][ipt]  = 0.3;   //ptK
        anacutsval[4][ipt]  = 0.3;   //ptPi
        anacutsval[5][ipt]  = 0.5;   //d0K
        anacutsval[6][ipt]  = 0.5;   //d0Pi
        anacutsval[7][ipt]  = 1.0;   //d0d0
        anacutsval[8][ipt]  = 0.7;   //cosThetaPoint
        anacutsval[9][ipt]  = 0.3;   //inv. mass half width of D* [GeV]
        anacutsval[10][ipt] = 0.3;   //half width of (M_Kpipi-M_D0) [GeV]
        anacutsval[11][ipt] = 0.05;  //PtMin of pi_s [GeV/c]
        anacutsval[13][ipt] = 0.5;   //theta, angle between the pi_s and decay plane of the D0 [rad]
    }

    //pT 3-5
    anacutsval[0][0] = 0.05;  //minv
    anacutsval[8][0] = 0.75;  //cosp
    anacutsval[12][0] = 1.;   //PtMax of pi_s
    anacutsval[14][0] = 0.75; //cosp
    anacutsval[15][0] = 1;    //NormDecayLenghtXY

    //pT 5-50
    anacutsval[0][1] = 0.10;  //minv
    anacutsval[8][1] = 0.80;  //cosp
    anacutsval[12][1] = 100.; //PtMax of pi_s
    anacutsval[14][1] = 0.80; //cospxy
    anacutsval[15][1] = 2;    //NormDecayLenghtXY


    AliRDHFCutsDStartoKpipi *analysiscuts = new AliRDHFCutsDStartoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dstar Analysis");
    analysiscuts->SetPtBins(nptbins + 1, ptbins);
    analysiscuts->SetCuts(nvars, nptbins, anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->AddTrackCutsSoftPi(esdTrackCutsSoftPi);

    //pid settings
    AliAODPidHF* pidObj = new AliAODPidHF();
    double priors[5] = {0.01, 0.001, 0.3, 0.3, 0.3};
    pidObj->SetPriors(priors, 5);
    pidObj->SetMatch(1);
    pidObj->SetSigma(0, 3); // TPC
    pidObj->SetSigma(3, 3); // TOF
    pidObj->SetTPC(true);
    pidObj->SetTOF(true);
    pidObj->SetOldPid(false);
    
    analysiscuts->SetPidHF(pidObj);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");
    if(triggername.EqualTo("kINT7"))
        analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    else if(triggername.EqualTo("kHighMultV0"))
        analysiscuts->SetTriggerMask(AliVEvent::kHighMultV0);
    else if(triggername.EqualTo("kINT7kHighMultV0"))
        analysiscuts->SetTriggerMask((AliVEvent::kHighMultV0 | AliVEvent::kINT7));
    else if(triggername.EqualTo("kMB"))
        analysiscuts->SetTriggerMask(AliVEvent::kMB);
    else {
        std::cout << "ERROR, trigger name " << triggername.Data() << " not supported! Exit" << std::endl;
        return nullptr;
    }

    analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
    analysiscuts->SetRemoveDaughtersFromPrim(true);
    analysiscuts->SetMinVtxContr(1);

    analysiscuts->SetMinPtCandidate(3.);
    analysiscuts->SetMaxPtCandidate(50.);

    std::cout << "This is the object I'm going to save:" << nptbins << std::endl;

    if (fIsMC)
        triggername = "kMB";

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = Form("DstartoKpipiCuts_pp13TeV_loose_%s.root", triggername.Data());
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}


AliRDHFCutsDStartoKpipi *MakeFileForCutsDstarpp13TeV_TreeML_Femto(bool fIsMC = false, TString triggername="kHighMultV0")
{

    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3, 1.e10);

    //soft pion selections
    AliESDtrackCuts* esdTrackCutsSoftPi = new AliESDtrackCuts();
    esdTrackCutsSoftPi->SetRequireSigmaToVertex(false);
    esdTrackCutsSoftPi->SetRequireTPCRefit(false);
    esdTrackCutsSoftPi->SetRequireITSRefit(false);
    esdTrackCutsSoftPi->SetMinNClustersITS(2);
    esdTrackCutsSoftPi->SetEtaRange(-0.8, +0.8);
    esdTrackCutsSoftPi->SetPtRange(0.05, 1.e10);

    const int nptbins = 2;
    float *ptbins;
    ptbins = new float[nptbins + 1];

    ptbins[0] = 1.;
    ptbins[1] = 5.;
    ptbins[2] = 50.;

    const int nvars = 16;
    float **anacutsval;
    anacutsval = new float *[nvars];
    for (int ic = 0; ic < nvars; ic++)
        anacutsval[ic] = new float[nptbins];

    /*
     Cut list
     0          "inv. mass [GeV]",
     1			"dca [cm]",
     2			"cosThetaStar",
     3			"pTK [GeV/c]",
     4			"pTPi [GeV/c]",
     5			"d0K [cm]",
     6			"d0Pi [cm]",
     7			"d0d0 [cm^2]",
     8			"cosThetaPoint",
     9			"inv. mass half width of D* [GeV]",
     10			"half width of (M_Kpipi-M_D0) [GeV]",
     11			"PtMin of pi_s [GeV/c]",
     12			"PtMax of pi_s [GeV/c]",
     13			"theta, angle between the pi_s and decay plane of the D0 [rad]",
     14			"|cosThetaPointXY|",
     15			"NormDecayLenghtXY";
     */

    for (int ipt = 0; ipt < nptbins; ipt++)
    {
        anacutsval[1][ipt]  = 0.1;   //dca
        anacutsval[2][ipt]  = 1.1;   //cost*
        anacutsval[3][ipt]  = 0.3;   //ptK
        anacutsval[4][ipt]  = 0.3;   //ptPi
        anacutsval[5][ipt]  = 0.5;   //d0K
        anacutsval[6][ipt]  = 0.5;   //d0Pi
        anacutsval[7][ipt]  = 1.0;   //d0d0
        anacutsval[8][ipt]  = 0.7;   //cosThetaPoint
        anacutsval[9][ipt]  = 0.3;   //inv. mass half width of D* [GeV]
        anacutsval[10][ipt] = 0.3;   //half width of (M_Kpipi-M_D0) [GeV]
        anacutsval[11][ipt] = 0.05;  //PtMin of pi_s [GeV/c]
        anacutsval[13][ipt] = 0.5;   //theta, angle between the pi_s and decay plane of the D0 [rad]
    }

    //pT 1-5
    anacutsval[0][0] = 0.05;  //minv
    anacutsval[8][0] = 0.75;  //cosp
    anacutsval[12][0] = 1.;   //PtMax of pi_s
    anacutsval[14][0] = 0.75; //cosp
    anacutsval[15][0] = 1;    //NormDecayLenghtXY

    //pT 5-50
    anacutsval[0][1] = 0.10;  //minv
    anacutsval[8][1] = 0.80;  //cosp
    anacutsval[12][1] = 100.; //PtMax of pi_s
    anacutsval[14][1] = 0.80; //cospxy
    anacutsval[15][1] = 2;    //NormDecayLenghtXY


    AliRDHFCutsDStartoKpipi *analysiscuts = new AliRDHFCutsDStartoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dstar Analysis");
    analysiscuts->SetPtBins(nptbins + 1, ptbins);
    analysiscuts->SetCuts(nvars, nptbins, anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->AddTrackCutsSoftPi(esdTrackCutsSoftPi);

    //pid settings
    AliAODPidHF* pidObj = new AliAODPidHF();
    double priors[5] = {0.01, 0.001, 0.3, 0.3, 0.3};
    pidObj->SetPriors(priors, 5);
    pidObj->SetMatch(1);
    pidObj->SetSigma(0, 3); // TPC
    pidObj->SetSigma(3, 3); // TOF
    pidObj->SetTPC(true);
    pidObj->SetTOF(true);
    pidObj->SetOldPid(false);
    
    analysiscuts->SetPidHF(pidObj);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");
    if(triggername.EqualTo("kINT7"))
        analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    else if(triggername.EqualTo("kHighMultV0"))
        analysiscuts->SetTriggerMask(AliVEvent::kHighMultV0);
    else if(triggername.EqualTo("kINT7kHighMultV0"))
        analysiscuts->SetTriggerMask((AliVEvent::kHighMultV0 | AliVEvent::kINT7));
    else if(triggername.EqualTo("kMB"))
        analysiscuts->SetTriggerMask(AliVEvent::kMB);
    else {
        std::cout << "ERROR, trigger name " << triggername.Data() << " not supported! Exit" << std::endl;
        return nullptr;
    }

    analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
    analysiscuts->SetRemoveDaughtersFromPrim(true);
    analysiscuts->SetMinVtxContr(1);

    analysiscuts->SetMinPtCandidate(1.);
    analysiscuts->SetMaxPtCandidate(50.);

    std::cout << "This is the object I'm going to save:" << nptbins << std::endl;

    if (fIsMC)
        triggername = "kMB";

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = Form("DstartoKpipiCuts_pp13TeV_loose_femto_%s.root", triggername.Data());
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

//__________________________________________________________________________________________
AliRDHFCutsDStartoKpipi *MakeFileForCutsDstarpp13TeV_NoDstarCuts(bool fIsMC = false, TString triggername="kINT7")
{

    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3, 1.e10);

    //soft pion selections
    AliESDtrackCuts* esdTrackCutsSoftPi = new AliESDtrackCuts();
    esdTrackCutsSoftPi->SetRequireSigmaToVertex(false);
    esdTrackCutsSoftPi->SetRequireTPCRefit(false);
    esdTrackCutsSoftPi->SetRequireITSRefit(false);
    esdTrackCutsSoftPi->SetMinNClustersITS(2);
    esdTrackCutsSoftPi->SetEtaRange(-0.8, +0.8);
    esdTrackCutsSoftPi->SetPtRange(0.05, 1.e10);

    const int nptbins = 2;
    float *ptbins;
    ptbins = new float[nptbins + 1];

    ptbins[0] = 1.;
    ptbins[1] = 5.;
    ptbins[2] = 50.;

    const int nvars = 16;
    float **anacutsval;
    anacutsval = new float *[nvars];
    for (int ic = 0; ic < nvars; ic++)
        anacutsval[ic] = new float[nptbins];

    /*
     Cut list
     0   "inv. mass [GeV]",
     1   "dca [cm]",
     2   "cosThetaStar",
     3   "pTK [GeV/c]",
     4   "pTPi [GeV/c]",
     5   "d0K [cm]",
     6   "d0Pi [cm]",
     7   "d0d0 [cm^2]",
     8   "cosThetaPoint",
     9   "inv. mass half width of D* [GeV]",
     10  "half width of (M_Kpipi-M_D0) [GeV]",
     11  "PtMin of pi_s [GeV/c]",
     12  "PtMax of pi_s [GeV/c]",
     13  "theta, angle between the pi_s and decay plane of the D0 [rad]",
     14  "|cosThetaPointXY|",
     15  "NormDecayLenghtXY";
    */

    for (int ipt = 0; ipt < nptbins; ipt++)
    {
        anacutsval[1][ipt]  = 0.1;   //dca O
        anacutsval[2][ipt]  = 1.1;   //cost* O
        anacutsval[3][ipt]  = 0.3;   //ptK  O
        anacutsval[4][ipt]  = 0.3;   //ptPi O
        anacutsval[5][ipt]  = 0.5;   //d0K  O
        anacutsval[6][ipt]  = 0.5;   //d0Pi O
        anacutsval[7][ipt]  = 1.0;   //d0d0 O
        anacutsval[9][ipt]  = 1000.; //inv. mass half width of D* [GeV]
        anacutsval[10][ipt] = 1000.; //half width of (M_Kpipi-M_D0) [GeV]
        anacutsval[11][ipt] = 0.05;  //PtMin of pi_s [GeV/c]
        anacutsval[12][ipt] = 1000.; //PtMax of pi_s [GeV/c]
        anacutsval[13][ipt] = 10.;   //theta, angle between the pi_s and decay plane of the D0 [rad]
    }

    //pT 1-5
    anacutsval[0][0] = 0.05;  //minv O
    anacutsval[8][0] = 0.75;  //cosp O
    anacutsval[14][0] = 0.75; //cosp xy O
    anacutsval[15][0] = 1;    //NormDecayLenghtXY O

    //pT 5-50
    anacutsval[0][1] = 0.10;  //minv
    anacutsval[8][1] = 0.80;  //cosp
    anacutsval[14][1] = 0.80; //cosp xy
    anacutsval[15][1] = 2;    //NormDecayLenghtXY


    AliRDHFCutsDStartoKpipi *analysiscuts = new AliRDHFCutsDStartoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dstar Analysis");
    analysiscuts->SetPtBins(nptbins + 1, ptbins);
    analysiscuts->SetCuts(nvars, nptbins, anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->AddTrackCutsSoftPi(esdTrackCutsSoftPi);

    //pid settings
    AliAODPidHF* pidObj = new AliAODPidHF();
    double priors[5] = {0.01, 0.001, 0.3, 0.3, 0.3};
    pidObj->SetPriors(priors, 5);
    pidObj->SetMatch(1);
    pidObj->SetSigma(0, 3); // TPC
    pidObj->SetSigma(3, 3); // TOF
    pidObj->SetTPC(true);
    pidObj->SetTOF(true);
    pidObj->SetOldPid(false);
    
    analysiscuts->SetPidHF(pidObj);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");
    if(triggername.EqualTo("kINT7"))
        analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    else if(triggername.EqualTo("kHighMultV0"))
        analysiscuts->SetTriggerMask(AliVEvent::kHighMultV0);
    else if(triggername.EqualTo("kINT7kHighMultV0"))
        analysiscuts->SetTriggerMask((AliVEvent::kHighMultV0 | AliVEvent::kINT7));
    else if(triggername.EqualTo("kMB"))
        analysiscuts->SetTriggerMask(AliVEvent::kMB);
    else {
        std::cout << "ERROR, trigger name " << triggername.Data() << " not supported! Exit" << std::endl;
        return nullptr;
    }

    analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
    analysiscuts->SetRemoveDaughtersFromPrim(true);
    analysiscuts->SetMinVtxContr(1);

    analysiscuts->SetMinPtCandidate(3.);
    analysiscuts->SetMaxPtCandidate(50.);

    std::cout << "This is the object I'm going to save:" << nptbins << std::endl;

    if (fIsMC)
        triggername = "kMB";

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = Form("DstartoKpipiCuts_pp13TeV_loose_NoDstarCuts_%s.root", triggername.Data());
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}
