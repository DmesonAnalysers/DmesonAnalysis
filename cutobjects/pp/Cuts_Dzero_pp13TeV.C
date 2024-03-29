#if !defined(__CINT__) || defined(__CLING__)

#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>

#include "AliESDtrackCuts.h"
#include "AliRDHFCutsD0toKpi.h"

#endif

//____________________________________________________________________________________________________//
// Methods:
// 1) MakeFileForCutsDzeropp13TeV --> loose cuts for 2023 Ds resonance analysis, matching Dstar ones
//____________________________________________________________________________________________________//

//__________________________________________________________________________________________
AliRDHFCutsD0toKpi *MakeFileForCutsDzeropp13TeV(bool fIsMC = false, TString triggername="kINT7")
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
    ptbins[1] = 5.;
    ptbins[2] = 50.;

    const int nvars = 11;
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
     9   "|cosThetaPointXY|",
     10  "NormDecayLenghtXY";
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
    }

    //pT 0-5
    anacutsval[0][0] = 0.05;  //minv
    anacutsval[8][0] = 0.75;  //cosp
    anacutsval[9][0] = 0.75;  //cosp xy
    anacutsval[10][0] = 1;    //NormDecayLenghtXY

    //pT 5-50
    anacutsval[0][1] = 0.10;  //minv
    anacutsval[8][1] = 0.80;  //cosp
    anacutsval[9][1] = 0.80;  //cosp xy
    anacutsval[10][1] = 2;    //NormDecayLenghtXY


    AliRDHFCutsD0toKpi *analysiscuts = new AliRDHFCutsD0toKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dzero Analysis");
    analysiscuts->SetPtBins(nptbins + 1, ptbins);
    analysiscuts->SetCuts(nvars, nptbins, anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);

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

    analysiscuts->SetMinPtCandidate(0.);
    analysiscuts->SetMaxPtCandidate(50.);

    std::cout << "This is the object I'm going to save:" << nptbins << std::endl;

    if (fIsMC)
        triggername = "kMB";

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = Form("DzerotoKpiCuts_pp13TeV_loose_%s.root", triggername.Data());
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}
