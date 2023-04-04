#if !defined(__CINT__) || defined(__CLING__)

#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>

#include "AliESDtrackCuts.h"
#include "AliRDHFCutsLctoV0.h"

#endif

//__________________________________________________________________________________________________________________//
// Methods:
// 1) MakeFileForCutsLc2V0bachelorpp5TeV_TreeML --> loose cuts for 2020 analysis (non-prompt Lc->pK0s and Lc->piL)
// 2) MakeFileForCutsLc2V0bachelorpp13TeV_TreeML --> loose cuts for 2020 analysis (non-prompt Lc->pK0s and Lc->piL)
//__________________________________________________________________________________________________________________//

//__________________________________________________________________________________________
AliRDHFCutsLctoV0 *MakeFileForCutsDspp13TeV_TreeML(bool fIsMC = false) {

  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(false);
  esdTrackCuts->SetRequireTPCRefit(true);
  esdTrackCuts->SetRequireITSRefit(true);
  esdTrackCuts->SetMinNCrossedRowsTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3, 1.e10);
  esdTrackCuts->SetEtaRange(-0.8, +0.8);
  esdTrackCuts->SetAcceptKinkDaughters(false);

  AliESDtrackCuts *esdTrackCutsV0daughters = new AliESDtrackCuts();
  esdTrackCutsV0daughters->SetRequireSigmaToVertex(false);
  esdTrackCutsV0daughters->SetRequireTPCRefit(true);
  esdTrackCutsV0daughters->SetRequireITSRefit(false);
  esdTrackCutsV0daughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCutsV0daughters->SetMinDCAToVertexXY(0.);
  esdTrackCutsV0daughters->SetPtRange(0.3, 1.e10);
  esdTrackCutsV0daughters->SetEtaRange(-0.8, +0.8);

  const int nptbins = 2;
  float *ptbins;
  ptbins = new float[nptbins + 1];

  ptbins[0] = 2.;
  ptbins[1] = 5.;
  ptbins[2] = 50.;

  const int nvars = 21;
  float **anacutsval;
  anacutsval = new float *[nvars];
  for (int ic = 0; ic < nvars; ic++)
    anacutsval[ic] = new float[nptbins];

  /*
     Cut list
        Lc inv. mass if K0S [GeV/c2]
        Lc inv. mass if Lambda [GeV/c2]
        K0S inv. mass [GeV/c2]
        Lambda/LambdaBar inv. mass[GeV/c2]
        pT min bachelor track [GeV/c]
        pT min V0-positive track [GeV/c]
        pT min V0-negative track [GeV/c]
        dca cascade (prong-to-prong) cut [cm]
        dca V0 (prong-to-prong) cut [number of sigmas]
        V0 cosPA min wrt PV
        d0 max bachelor wrt PV [cm]
        d0 max V0 wrt PV [cm]
        mass K0S veto [GeV/c2]
        mass Lambda/LambdaBar veto [GeV/c2]
        mass Gamma veto [GeV/c2]
        pT min V0 track [GeV/c]
        Max Proton emission angle in Lc CMS
        Min Proton emission angle in Lc CMS
        Resigned d0
        V0 qT/|alpha|
        V0 type
    */

  for (int ipt = 0; ipt < nptbins; ipt++) {
    anacutsval[0][ipt] = 0.2;
    anacutsval[1][ipt] = 0.2;
    anacutsval[2][ipt] = 0.03;
    anacutsval[3][ipt] = 0.05;
    anacutsval[4][ipt] = 0.0;
    anacutsval[5][ipt] = 0.0;
    anacutsval[6][ipt] = 0.0;
    anacutsval[7][ipt] = 1000.;
    anacutsval[8][ipt] = 0.8;
    anacutsval[9][ipt] = 0.997;
    anacutsval[10][ipt] = 3.0;
    anacutsval[11][ipt] = 1000;
    anacutsval[12][ipt] = 0.0;
    anacutsval[13][ipt] = 0.0;
    anacutsval[14][ipt] = 0.0;
    anacutsval[15][ipt] = 0.0;
    anacutsval[16][ipt] = 0.9;
    anacutsval[17][ipt] = -0.9;
    anacutsval[18][ipt] = -0.8;
    anacutsval[19][ipt] = 1.8;
    anacutsval[20][ipt] = 0;
  }

  AliRDHFCutsLctoV0 *analysiscuts = new AliRDHFCutsLctoV0();
  analysiscuts->SetName("AnalysisCuts");

  analysiscuts->SetTitle("Cuts for non-prompt Lc->V0bachelor analysis");
  analysiscuts->SetPtBins(nptbins + 1, ptbins);
  analysiscuts->SetCuts(nvars, nptbins, anacutsval);
  analysiscuts->AddTrackCuts(esdTrackCuts);
  analysiscuts->SetKinkRejection(!esdTrackCuts->GetAcceptKinkDaughters());
  analysiscuts->AddTrackCutsV0daughters(esdTrackCutsV0daughters);
  analysiscuts->SetUseTrackSelectionWithFilterBits(false);

  // PID preselection for bachelor (3sigma)
  AliAODPidHF *pidObjBachelor = new AliAODPidHF();
  double sigmasBac[5] = {3., 3., 3., 3., 3.}; // 0, 1(A), 2(A) -> TPC; 3 -> TOF; 4 -> ITS
  pidObjBachelor->SetSigma(sigmasBac);
  pidObjBachelor->SetAsym(false);
  pidObjBachelor->SetMatch(1);
  pidObjBachelor->SetTPC(true);
  pidObjBachelor->SetTOF(true);
  pidObjBachelor->SetTOFdecide(false);

  analysiscuts->SetPidSelectionFlag(11);
  analysiscuts->SetPidHF(pidObjBachelor);
  analysiscuts->SetUsePID(true);

  analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
  analysiscuts->SetTriggerClass("");
  analysiscuts->SetTriggerMask(AliVEvent::kINT7);
  if (fIsMC)
    analysiscuts->SetTriggerMask(AliVEvent::kMB);

  analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
  analysiscuts->SetRemoveDaughtersFromPrim(true);
  analysiscuts->SetMinVtxContr(1);

  analysiscuts->SetMinPtCandidate(2.);
  analysiscuts->SetMaxPtCandidate(50.);

  TString triggername = "kINT7";
  if (fIsMC)
    triggername = "kMB";

  analysiscuts->PrintAll();
  analysiscuts->PrintTrigger();
  TString filename = Form("LctoV0bachelorCuts_pp_13TeV_nonprompt_loose_%s.root", triggername.Data());
  TFile *fout = new TFile(filename.Data(), "recreate");
  fout->cd();
  analysiscuts->Write();
  fout->Close();

  return analysiscuts;
}


//__________________________________________________________________________________________
AliRDHFCutsLctoV0 *MakeFileForCutsDspp5TeV_TreeML(bool fIsMC = false)
{

    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3, 1.e10);
    esdTrackCuts->SetEtaRange(-0.8,+0.8);
    esdTrackCuts->SetAcceptKinkDaughters(false);

    AliESDtrackCuts *esdTrackCutsV0daughters = new AliESDtrackCuts();
    esdTrackCutsV0daughters->SetRequireSigmaToVertex(false);
    esdTrackCutsV0daughters->SetRequireTPCRefit(true);
    esdTrackCutsV0daughters->SetRequireITSRefit(false);
    esdTrackCutsV0daughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCutsV0daughters->SetMinDCAToVertexXY(0.);
    esdTrackCutsV0daughters->SetPtRange(0.3, 1.e10);
    esdTrackCutsV0daughters->SetEtaRange(-0.8,+0.8);

    const int nptbins = 2;
    float *ptbins;
    ptbins = new float[nptbins + 1];

    ptbins[0] = 2.;
    ptbins[1] = 5.;
    ptbins[2] = 50.;

    const int nvars = 21;
    float **anacutsval;
    anacutsval = new float *[nvars];
    for (int ic = 0; ic < nvars; ic++)
        anacutsval[ic] = new float[nptbins];

    /*
     Cut list
        Lc inv. mass if K0S [GeV/c2]
        Lc inv. mass if Lambda [GeV/c2]
        K0S inv. mass [GeV/c2]
        Lambda/LambdaBar inv. mass[GeV/c2]
        pT min bachelor track [GeV/c]
        pT min V0-positive track [GeV/c]
        pT min V0-negative track [GeV/c]
        dca cascade (prong-to-prong) cut [cm]
        dca V0 (prong-to-prong) cut [number of sigmas]
        V0 cosPA min wrt PV
        d0 max bachelor wrt PV [cm]
        d0 max V0 wrt PV [cm]
        mass K0S veto [GeV/c2]
        mass Lambda/LambdaBar veto [GeV/c2]
        mass Gamma veto [GeV/c2]
        pT min V0 track [GeV/c]
        Max Proton emission angle in Lc CMS
        Min Proton emission angle in Lc CMS
        Resigned d0
        V0 qT/|alpha|
        V0 type
    */

    for (int ipt = 0; ipt < nptbins; ipt++)
    {
        anacutsval[0][ipt]=0.2;
        anacutsval[1][ipt]=0.2;
        anacutsval[2][ipt]=0.03;
        anacutsval[3][ipt]=0.05;
        anacutsval[4][ipt]=0.0;
        anacutsval[5][ipt]=0.0;
        anacutsval[6][ipt]=0.0;
        anacutsval[7][ipt]=1000.;
        anacutsval[8][ipt]=0.8;
        anacutsval[9][ipt]=0.997;
        anacutsval[10][ipt]=3.0;
        anacutsval[11][ipt]=1.5;
        anacutsval[12][ipt]=0.0;
        anacutsval[13][ipt]=0.0;
        anacutsval[14][ipt]=0.0;
        anacutsval[15][ipt]=0.0;
        anacutsval[16][ipt]=0.9;
        anacutsval[17][ipt]=-0.9;
        anacutsval[18][ipt]=-0.8;
        anacutsval[19][ipt]=1.8;
        anacutsval[20][ipt]=0;
    }

    AliRDHFCutsLctoV0* analysiscuts = new AliRDHFCutsLctoV0();
    analysiscuts->SetName("AnalysisCuts");
    
    analysiscuts->SetTitle("Cuts for non-prompt Lc->V0bachelor analysis");
    analysiscuts->SetPtBins(nptbins + 1, ptbins);
    analysiscuts->SetCuts(nvars, nptbins, anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetKinkRejection(!esdTrackCuts->GetAcceptKinkDaughters());
    analysiscuts->AddTrackCutsV0daughters(esdTrackCutsV0daughters);
    analysiscuts->SetUseTrackSelectionWithFilterBits(false);
    
    // PID preselection for bachelor (3sigma)    
    AliAODPidHF* pidObjBachelor = new AliAODPidHF();
    double sigmasBac[5] = {3., 3., 3., 3., 3.}; // 0, 1(A), 2(A) -> TPC; 3 -> TOF; 4 -> ITS
    pidObjBachelor->SetSigma(sigmasBac);
    pidObjBachelor->SetAsym(false);
    pidObjBachelor->SetMatch(1);
    pidObjBachelor->SetTPC(true);
    pidObjBachelor->SetTOF(true);
    pidObjBachelor->SetTOFdecide(false);
    
    analysiscuts->SetPidSelectionFlag(11);
    analysiscuts->SetPidHF(pidObjBachelor);
    analysiscuts->SetUsePID(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");
    analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    if (fIsMC)
        analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
    analysiscuts->SetRemoveDaughtersFromPrim(true);
    analysiscuts->SetMinVtxContr(1);

    analysiscuts->SetMinPtCandidate(2.);
    analysiscuts->SetMaxPtCandidate(50.);

    TString triggername = "kINT7";
    if (fIsMC)
        triggername = "kMB";

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = Form("LctoV0bachelorCuts_pp_nonprompt_loose_%s.root", triggername.Data());
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}
