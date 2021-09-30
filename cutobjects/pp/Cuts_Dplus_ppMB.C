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
// 1) MakeFileForCutsDpluspp5TeV_CompO2 --> selections from 2017 analysis 
//                                          https://alice-notes.web.cern.ch/node/808 adapted for O2 comparison
//                                          
//____________________________________________________________________________________________________//

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi *MakeFileForCutsDpluspp5TeV_CompO2(bool fIsMC = false)
{

    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3, 1.e10);

    const int nptbins = 12;
    float *ptbins;
    ptbins = new float[nptbins + 1];

    ptbins[0]=1.;
    ptbins[1]=2.;
    ptbins[2]=3.;
    ptbins[3]=4.;
    ptbins[4]=5.;
    ptbins[5]=6.;
    ptbins[6]=7.;
    ptbins[7]=8.;
    ptbins[8]=10.;
    ptbins[9]=12.;
    ptbins[10]=16.;
    ptbins[11]=24.;
    ptbins[12]=36.;

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

    for (int ipt = 0; ipt < nptbins; ipt++){
        anacutsval[0][ipt]=0.2;    //minv
        anacutsval[1][ipt]=0.3;    //ptK
        anacutsval[2][ipt]=0.3;    //ptPi
        anacutsval[3][ipt]=0.;     //d0K
        anacutsval[4][ipt]=0.;     //d0Pi
        anacutsval[5][ipt]=0.;     //dist12
        anacutsval[6][ipt]=1000.;  //sigvert
        anacutsval[8][ipt]=0.0;    //ptMax
        anacutsval[10][ipt]=0.0;   //sumd02
        anacutsval[11][ipt]=1000.; //dca
    }

    //declen
    anacutsval[7][0]=0.07;   //1.0-2.0
    anacutsval[7][1]=0.07;   //2.0-3.0
    anacutsval[7][2]=0.10;   //3.0-4.0
    anacutsval[7][3]=0.10;   //4.0-5.0
    anacutsval[7][4]=0.10;   //5.0-6.0
    anacutsval[7][5]=0.10;   //6.0-7.0
    anacutsval[7][6]=0.10;   //7.0-8.0
    anacutsval[7][7]=0.12;   //8.0-10.0
    anacutsval[7][8]=0.12;   //10.0-12.0
    anacutsval[7][9]=0.12;   //12.0-16.0
    anacutsval[7][10]=0.12;  //16.0-24.0
    anacutsval[7][11]=0.20;  //24.0-36.0


    //cosp
    anacutsval[9][0]=0.96;   //1.0-2.0
    anacutsval[9][1]=0.96;   //2.0-3.0
    anacutsval[9][2]=0.96;   //3.0-4.0
    anacutsval[9][3]=0.96;   //4.0-5.0
    anacutsval[9][4]=0.96;   //5.0-6.0
    anacutsval[9][5]=0.96;   //6.0-7.0
    anacutsval[9][6]=0.96;   //7.0-8.0
    anacutsval[9][7]=0.96;   //8.0-10.0
    anacutsval[9][8]=0.96;   //10.0-12.0
    anacutsval[9][9]=0.96;   //12.0-16.0
    anacutsval[9][10]=0.96;  //16.0-24.0
    anacutsval[9][11]=0.94;  //24.0-36.0

    //ndlXY
    anacutsval[12][0]=6.;   //1.0-2.0
    anacutsval[12][1]=5.;   //2.0-3.0
    anacutsval[12][2]=5.;   //3.0-4.0
    anacutsval[12][3]=5.;   //4.0-5.0
    anacutsval[12][4]=5.;   //5.0-6.0
    anacutsval[12][5]=5.;   //6.0-7.0
    anacutsval[12][6]=5.;   //7.0-8.0
    anacutsval[12][7]=5.;   //8.0-10.0
    anacutsval[12][8]=5.;   //10.0-12.0
    anacutsval[12][9]=5.;   //12.0-16.0
    anacutsval[12][10]=5.;  //16.0-24.0
    anacutsval[12][11]=5.;  //24.0-36.0

    //cospXY
    anacutsval[13][0]=0.985;   //1.0-2.0
    anacutsval[13][1]=0.985;   //2.0-3.0
    anacutsval[13][2]=0.980;   //3.0-4.0
    anacutsval[13][3]=0.000;   //4.0-5.0
    anacutsval[13][4]=0.000;   //5.0-6.0
    anacutsval[13][5]=0.000;   //6.0-7.0
    anacutsval[13][6]=0.000;   //7.0-8.0
    anacutsval[13][7]=0.000;   //8.0-10.0
    anacutsval[13][8]=0.000;   //10.0-12.0
    anacutsval[13][9]=0.000;   //12.0-16.0
    anacutsval[13][10]=0.000;  //16.0-24.0
    anacutsval[13][11]=0.000;  //24.0-36.0

    //d0d0exp
    float *d0d0expcutsval = new float[nptbins];
    for (int ipt = 0; ipt < nptbins; ipt++){
        d0d0expcutsval[ipt]=2.5;
    }

    AliRDHFCutsDplustoKpipi *analysiscuts = new AliRDHFCutsDplustoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
    analysiscuts->SetPtBins(nptbins + 1, ptbins);
    analysiscuts->SetCuts(nvars, nptbins, anacutsval);
    analysiscuts->Setd0MeasMinusExpCut(nptbins, d0d0expcutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetScaleNormDLxyBypOverPt(false);

    analysiscuts->SetUsePID(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");
    analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    if (fIsMC)
        analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUseImpParProdCorrCut(false);

    analysiscuts->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
    analysiscuts->SetRemoveDaughtersFromPrim(true);
    analysiscuts->SetMinVtxContr(1);

    analysiscuts->SetMinPtCandidate(1.);
    analysiscuts->SetMaxPtCandidate(36.);

    std::cout << "This is the object I'm going to save:" << nptbins << std::endl;

    TString triggername = "kINT7";
    if (fIsMC)
        triggername = "kMB";

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename = Form("DplustoKpipiCuts_pp_compO2_%s.root", triggername.Data());
    TFile *fout = new TFile(filename.Data(), "recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}
