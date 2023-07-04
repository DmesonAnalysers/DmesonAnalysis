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
// 1) MakeFileForCutsDpluspPb5TeV_TreeML --> loose cuts for 2023 analysis
//____________________________________________________________________________________________________//

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi *MakeFileForCutsDpluspPb5TeV_TreeML()
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

    ptbins[0] = 2.;
    ptbins[1] = 5.;
    ptbins[2] = 50.;

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
    }

    //pT 2-5
    anacutsval[6][0] = 0.040; //sigvert
    anacutsval[7][0] = 0.030; //declen
    anacutsval[9][0] = 0.90;  //cosp
    anacutsval[12][0] = 4.;   //ndlXY
    anacutsval[13][0] = 0.90; //cospXY

    //pT 5-50
    anacutsval[6][1] = 0.050; //sigvert
    anacutsval[7][1] = 0.040; //declen
    anacutsval[9][1] = 0.85;  //cosp
    anacutsval[12][1] = 3.;   //ndlXY
    anacutsval[13][1] = 0.85; //cospXY

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
    analysiscuts->SetTriggerMask(AliVEvent::kINT7); // high multiplicity V0M
    analysiscuts->SetUseImpParProdCorrCut(false);

    analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
    analysiscuts->SetRemoveDaughtersFromPrim(true);

    analysiscuts->SetMinPtCandidate(2.);
    analysiscuts->SetMaxPtCandidate(50.);

    std::cout << "This is the object I'm going to save:" << nptbins << std::endl;

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = "DplustoKpipiCuts_pPb_loose_kINT7.root";
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}
