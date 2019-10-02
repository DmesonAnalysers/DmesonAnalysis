#if !defined (__CINT__) || defined (__CLING__)

#include <iostream>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>

#include "AliESDtrackCuts.h"
#include "AliRDHFCutsDstoKKpi.h"

#endif

//____________________________________________________________________________________________________//
// Methods:
// 1) MakeFileForCutsDs6080_Central2018 --> central cuts of 2018 analysis
// 2) MakeFileForCutsDs6080_Loose2018 --> loose cuts of 2018 analysis (sparse)
// 3) MakeFileForCutsDs6080_FiltTreeCreator2018 --> filtering cuts of 2018 analysis (tree creator)
//____________________________________________________________________________________________________//

AliRDHFCutsDstoKKpi* MakeFileForCutsDs6080_Central2018(bool fUseStrongPID = true, double maxPtstrongPID = 8.0, bool fIsMC=false, int addTrackCut = 0, double ptmin = 2., double ptmax = 50.) {

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3,1.e10);

    float mincen=60.;
    float maxcen=80.;

    if(fIsMC){
      mincen=0.;
      maxcen=100.;
    }

    const int nptbins=9;
    float* ptbins;
    ptbins=new float[nptbins+1];
    ptbins[0]=2.;
    ptbins[1]=3.;
    ptbins[2]=4.;
    ptbins[3]=5.;
    ptbins[4]=6.;
    ptbins[5]=8.;
    ptbins[6]=12.;
    ptbins[7]=16.;
    ptbins[8]=24.;
    ptbins[9]=36.;

    const int nvars=20;

    float** anacutsval;
    anacutsval=new float*[nvars];

    for(int ic=0;ic<nvars;ic++){anacutsval[ic]=new float[nptbins];}
    for(int ipt=0;ipt<nptbins;ipt++){

        anacutsval[0][ipt]=0.2;
        anacutsval[1][ipt]=0.3;
        anacutsval[2][ipt]=0.3;
        anacutsval[3][ipt]=0.;
        anacutsval[4][ipt]=0.;
        anacutsval[5][ipt]=0.005;
        anacutsval[8][ipt]=0.;
        anacutsval[10][ipt]=0.;
        anacutsval[11][ipt]=1000.0;
        anacutsval[13][ipt]=0.001;
    }
    /*

     Cut list                                           rejection condition
     0      "inv. mass [GeV]",                          invmassDS-massDspdg>fCutsRD
     1			"pTK [GeV/c]",                              pTK<fCutsRd
     2			"pTPi [GeV/c]",                             pTPi<fCutsRd
     3			"d0K [cm]",                                 d0K<fCutsRd
     4			"d0Pi [cm]",                                d0Pi<fCutsRd
     5			"dist12 [cm]",                              dist12<fCutsRd
     6			"sigmavert [cm]",                           sigmavert>fCutsRd
     7			"decLen [cm]",                              decLen<fCutsRD
     8			"ptMax [GeV/c]",                            ptMax<fCutsRD
     9			"cosThetaPoint",                            CosThetaPoint<fCutsRD
     10			"Sum d0^2 (cm^2)",                          sumd0<fCutsRD
     11			"dca [cm]",                                 dca(i)>fCutsRD
     12			"inv. mass (Mphi-MKK) [GeV]",               invmass-pdg>fCutsRD
     13			"inv. mass (MKo*-MKpi) [GeV]"};             invmass-pdg>fCutsRD
     14    	"Abs(CosineKpiPhiRFrame)^3",
     15  		"CosPiDsLabFrame"};
     16  		"DecLengthXY
     17  		"NormDecayLength"};
     18  		"NormDecayLengthXY"};
     19  		"cosThetaPointXY"};
     */

    anacutsval[6][0]=0.030;   //sigmavert
    anacutsval[7][0]=0.03;   //decay length
    anacutsval[9][0]=0.99;   //cosP
    anacutsval[12][0]=0.008; //Mass Phi
    anacutsval[14][0]=0.05;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][0]=1.0;  //CosP labFrame
    anacutsval[16][0]=0.;  //decayXY
    anacutsval[17][0]=0.;    //normdecay
    anacutsval[18][0]=7.0;   //normdecayXY
    anacutsval[19][0]=0.99;  //CosPXY

    anacutsval[6][1]=0.035;   //sigmavert
    anacutsval[7][1]=0.03;   //decay length
    anacutsval[9][1]=0.97;   //cosP
    anacutsval[12][1]=0.008; //Mass Phi
    anacutsval[14][1]=0.0;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][1]=1.0;  //CosP labFrame
    anacutsval[16][1]=0.;  //decayXY
    anacutsval[17][1]=0.;    //normdecay
    anacutsval[18][1]=5.0;   //normdecayXY
    anacutsval[19][1]=0.97;  //CosPXY

    anacutsval[6][2]=0.035;   //sigmavert
    anacutsval[7][2]=0.03;   //decay length
    anacutsval[9][2]=0.97;   //cosP
    anacutsval[12][2]=0.008; //Mass Phi
    anacutsval[14][2]=0.0;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][2]=1.0;   //CosP labFrame
    anacutsval[16][2]=0.;  //decayXY
    anacutsval[17][2]=0.;    //normdecay
    anacutsval[18][2]=4.0;   //normdecayXY
    anacutsval[19][2]=0.97;  //CosPXY

    anacutsval[6][3]=0.035;   //sigmavert
    anacutsval[7][3]=0.03;   //decay length
    anacutsval[9][3]=0.97;   //cosP
    anacutsval[12][3]=0.008; //Mass Phi
    anacutsval[14][3]=0.0;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][3]=1.0;  //CosP labFrame
    anacutsval[16][3]=0.;   //decayXY
    anacutsval[17][3]=0.;    //normdecay
    anacutsval[18][3]=4.0;   //normdecayXY
    anacutsval[19][3]=0.97;  //CosPXY

    anacutsval[6][4]=0.035;   //sigmavert
    anacutsval[7][4]=0.03;   //decay length
    anacutsval[9][4]=0.97;   //cosP
    anacutsval[12][4]=0.008; //Mass Phi
    anacutsval[14][4]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][4]=1.0;  //CosP labFrame
    anacutsval[16][4]=0.;   //decayXY
    anacutsval[17][4]=0.;    //normdecay
    anacutsval[18][4]=4.0;   //normdecayXY
    anacutsval[19][4]=0.97;  //CosPXY

    anacutsval[6][5]=0.035;   //sigmavert
    anacutsval[7][5]=0.03;   //decay length
    anacutsval[9][5]=0.97;   //cosP
    anacutsval[12][5]=0.008; //Mass Phi
    anacutsval[14][5]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][5]=1.0;  //CosP labFrame
    anacutsval[16][5]=0.;   //decayXY
    anacutsval[17][5]=0.;    //normdecay
    anacutsval[18][5]=4.0;   //normdecayXY
    anacutsval[19][5]=0.97;  //CosPXY

    anacutsval[6][6]=0.050;   //sigmavert
    anacutsval[7][6]=0.03;   //decay length
    anacutsval[9][6]=0.98;   //cosP
    anacutsval[12][6]=0.010; //Mass Phi
    anacutsval[14][6]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][6]=1.0;  //CosP labFrame
    anacutsval[16][6]=0.;   //decayXY
    anacutsval[17][6]=0.;    //normdecay
    anacutsval[18][6]=4.0;   //normdecayXY
    anacutsval[19][6]=0.98;  //CosPXY

    anacutsval[6][7]=0.050;   //sigmavert
    anacutsval[7][7]=0.03;   //decay length
    anacutsval[9][7]=0.99;   //cosP
    anacutsval[12][7]=0.008; //Mass Phi
    anacutsval[14][7]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][7]=1.0;  //CosP labFrame
    anacutsval[16][7]=0.;   //decayXY
    anacutsval[17][7]=0.;    //normdecay
    anacutsval[18][7]=4.0;   //normdecayXY
    anacutsval[19][7]=0.99;  //CosPXY

    anacutsval[6][8]=0.050;   //sigmavert
    anacutsval[7][8]=0.03;   //decay length
    anacutsval[9][8]=0.99;   //cosP
    anacutsval[12][8]=0.008; //Mass Phi
    anacutsval[14][8]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][8]=1.0;  //CosP labFrame
    anacutsval[16][8]=0.;   //decayXY
    anacutsval[17][8]=0.;    //normdecay
    anacutsval[18][8]=4.0;   //normdecayXY
    anacutsval[19][8]=0.99;  //CosPXY

    float topomCuts[nptbins] = {3.,3.,3.,3.,3.,3.,3.,3.,3.}; //topomatic
    float d0Cuts[nptbins] = {1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.,1000.}; //impparXY

    AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Ds Analysis and CF");
    analysiscuts->SetUsePreSelect(1);
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->Setd0MeasMinusExpCut(nptbins,topomCuts);
    analysiscuts->Setd0Cut(nptbins,d0Cuts);
    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    if(!fIsMC)
      analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    else
      analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUsePID(true);
    if(fUseStrongPID) {
      analysiscuts->SetPidOption(1); //0=kConservative,1=kStrong
      analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
    }
    else analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb6080);

    analysiscuts->SetOptPileup(false);
    analysiscuts->SetMinCentrality(mincen);
    analysiscuts->SetMaxCentrality(maxcen);

    analysiscuts->SetMinPtCandidate(ptmin);
    analysiscuts->SetMaxPtCandidate(ptmax);
    analysiscuts->SetRemoveDaughtersFromPrim(false);

    TString trackCutName= "";
    switch (addTrackCut) {
      case 1:
        esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
        trackCutName = "_addRowsOverClusterTPC";
        break;

      case 2:
        analysiscuts->SetMinCrossedRowsTPCPtDep("120-(5/pt)");
        trackCutName = "_addMinCrossedRowsTPC";
        break;

      case 3:
        analysiscuts->SetMinRatioClsOverCrossRowsTPC(0.65);
        trackCutName = "_addRatioClsOverRowsTPC";
        break;

      case 4:
        analysiscuts->SetUseCutGeoNcrNcl(true);
        trackCutName = "_addUseCutGeo";
        break;

      default:
        break;
    }

    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;

    analysiscuts->PrintAll();
    TString triggername = "kINT7";
    if(fIsMC) triggername = "kMB";
    TString pidname = "";
    if(fUseStrongPID) pidname = Form("_strongPIDpt%0.f", maxPtstrongPID);

    TFile* fout=new TFile(Form("DstoKKpiCuts_6080_central%s_Raa_%s%s_pt%0.f_%0.f.root", pidname.Data(), triggername.Data(), trackCutName.Data(), ptmin, ptmax),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

AliRDHFCutsDstoKKpi* MakeFileForCutsDs6080_Loose2018(bool fUseStrongPID = true, double maxPtstrongPID = 8.0, bool fIsMC=false, double ptmin=2., double ptmax=50.) {

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3,1.e10);

    float mincen=60.;
    float maxcen=80.;

    if(fIsMC){
      mincen=0.;
      maxcen=100.;
    }

    const int nptbins=9;
    float* ptbins;
    ptbins=new float[nptbins+1];
    ptbins[0]=2.;
    ptbins[1]=3.;
    ptbins[2]=4.;
    ptbins[3]=5.;
    ptbins[4]=6.;
    ptbins[5]=8.;
    ptbins[6]=12.;
    ptbins[7]=16.;
    ptbins[8]=24.;
    ptbins[9]=36.;

    const int nvars=20;

    float** anacutsval;
    anacutsval=new float*[nvars];

    for(int ic=0;ic<nvars;ic++){anacutsval[ic]=new float[nptbins];}
    for(int ipt=0;ipt<nptbins;ipt++){

        anacutsval[0][ipt]=0.2;
        anacutsval[1][ipt]=0.3;
        anacutsval[2][ipt]=0.3;
        anacutsval[3][ipt]=0.;
        anacutsval[4][ipt]=0.;
        anacutsval[5][ipt]=0.005;
        anacutsval[8][ipt]=0.;
        anacutsval[10][ipt]=0.;
        anacutsval[11][ipt]=1000.0;
        anacutsval[13][ipt]=0.001;
    }
    /*

     Cut list                                           rejection condition
     0      "inv. mass [GeV]",                          invmassDS-massDspdg>fCutsRD
     1			"pTK [GeV/c]",                              pTK<fCutsRd
     2			"pTPi [GeV/c]",                             pTPi<fCutsRd
     3			"d0K [cm]",                                 d0K<fCutsRd
     4			"d0Pi [cm]",                                d0Pi<fCutsRd
     5			"dist12 [cm]",                              dist12<fCutsRd
     6			"sigmavert [cm]",                           sigmavert>fCutsRd
     7			"decLen [cm]",                              decLen<fCutsRD
     8			"ptMax [GeV/c]",                            ptMax<fCutsRD
     9			"cosThetaPoint",                            CosThetaPoint<fCutsRD
     10			"Sum d0^2 (cm^2)",                          sumd0<fCutsRD
     11			"dca [cm]",                                 dca(i)>fCutsRD
     12			"inv. mass (Mphi-MKK) [GeV]",               invmass-pdg>fCutsRD
     13			"inv. mass (MKo*-MKpi) [GeV]"};             invmass-pdg>fCutsRD
     14    	"Abs(CosineKpiPhiRFrame)^3",
     15  		"CosPiDsLabFrame"};
     16  		"DecLengthXY
     17  		"NormDecayLength"};
     18  		"NormDecayLengthXY"};
     19  		"cosThetaPointXY"};
     */

    anacutsval[6][0]=0.045;   //sigmavert
    anacutsval[7][0]=0.02;   //decay length
    anacutsval[9][0]=0.96;   //cosP
    anacutsval[12][0]=0.010; //Mass Phi
    anacutsval[14][0]=0.0;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][0]=1.0;  //CosP labFrame
    anacutsval[16][0]=0.;  //decayXY
    anacutsval[17][0]=0.;    //normdecay
    anacutsval[18][0]=3.0;   //normdecayXY
    anacutsval[19][0]=0.96;  //CosPXY

    anacutsval[6][1]=0.045;   //sigmavert
    anacutsval[7][1]=0.02;   //decay length
    anacutsval[9][1]=0.96;   //cosP
    anacutsval[12][1]=0.010; //Mass Phi
    anacutsval[14][1]=0.0;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][1]=1.0;  //CosP labFrame
    anacutsval[16][1]=0.;  //decayXY
    anacutsval[17][1]=0.;    //normdecay
    anacutsval[18][1]=3.0;   //normdecayXY
    anacutsval[19][1]=0.96;  //CosPXY

    anacutsval[6][2]=0.045;   //sigmavert
    anacutsval[7][2]=0.02;   //decay length
    anacutsval[9][2]=0.96;   //cosP
    anacutsval[12][2]=0.008; //Mass Phi
    anacutsval[14][2]=0.0;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][2]=1.0;   //CosP labFrame
    anacutsval[16][2]=0.;  //decayXY
    anacutsval[17][2]=0.;    //normdecay
    anacutsval[18][2]=3.0;   //normdecayXY
    anacutsval[19][2]=0.96;  //CosPXY

    anacutsval[6][3]=0.045;   //sigmavert
    anacutsval[7][3]=0.02;   //decay length
    anacutsval[9][3]=0.96;   //cosP
    anacutsval[12][3]=0.010; //Mass Phi
    anacutsval[14][3]=0.0;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][3]=1.0;  //CosP labFrame
    anacutsval[16][3]=0.;   //decayXY
    anacutsval[17][3]=0.;    //normdecay
    anacutsval[18][3]=3.0;   //normdecayXY
    anacutsval[19][3]=0.96;  //CosPXY

    anacutsval[6][4]=0.045;   //sigmavert
    anacutsval[7][4]=0.02;   //decay length
    anacutsval[9][4]=0.95;   //cosP
    anacutsval[12][4]=0.010; //Mass Phi
    anacutsval[14][4]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][4]=1.0;  //CosP labFrame
    anacutsval[16][4]=0.;   //decayXY
    anacutsval[17][4]=0.;    //normdecay
    anacutsval[18][4]=3.0;   //normdecayXY
    anacutsval[19][4]=0.95;  //CosPXY

    anacutsval[6][5]=0.045;   //sigmavert
    anacutsval[7][5]=0.02;   //decay length
    anacutsval[9][5]=0.95;   //cosP
    anacutsval[12][5]=0.010; //Mass Phi
    anacutsval[14][5]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][5]=1.0;  //CosP labFrame
    anacutsval[16][5]=0.;   //decayXY
    anacutsval[17][5]=0.;    //normdecay
    anacutsval[18][5]=3.0;   //normdecayXY
    anacutsval[19][5]=0.95;  //CosPXY

    anacutsval[6][6]=0.055;   //sigmavert
    anacutsval[7][6]=0.02;   //decay length
    anacutsval[9][6]=0.95;   //cosP
    anacutsval[12][6]=0.010; //Mass Phi
    anacutsval[14][6]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][6]=1.0;  //CosP labFrame
    anacutsval[16][6]=0.;   //decayXY
    anacutsval[17][6]=0.;    //normdecay
    anacutsval[18][6]=3.0;   //normdecayXY
    anacutsval[19][6]=0.95;  //CosPXY

    anacutsval[6][7]=0.055;   //sigmavert
    anacutsval[7][7]=0.02;   //decay length
    anacutsval[9][7]=0.95;   //cosP
    anacutsval[12][7]=0.010; //Mass Phi
    anacutsval[14][7]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][7]=1.0;  //CosP labFrame
    anacutsval[16][7]=0.;   //decayXY
    anacutsval[17][7]=0.;    //normdecay
    anacutsval[18][7]=3.0;   //normdecayXY
    anacutsval[19][7]=0.95;  //CosPXY

    anacutsval[6][8]=0.055;   //sigmavert
    anacutsval[7][8]=0.02;   //decay length
    anacutsval[9][8]=0.95;   //cosP
    anacutsval[12][8]=0.010; //Mass Phi
    anacutsval[14][8]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][8]=1.0;  //CosP labFrame
    anacutsval[16][8]=0.;   //decayXY
    anacutsval[17][8]=0.;    //normdecay
    anacutsval[18][8]=3.0;   //normdecayXY
    anacutsval[19][8]=0.95;  //CosPXY

    float topomCuts[nptbins] = {4.,4.,4.,5.,5.,5.,6.,6.,6.}; //topomatic
    float d0Cuts[nptbins] = {500.,500.,500.,500.,500.,500.,500.,500.,500.}; //impparXY

    AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Ds Analysis and CF");
    analysiscuts->SetUsePreSelect(1);
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->Setd0MeasMinusExpCut(nptbins,topomCuts);
    analysiscuts->Setd0Cut(nptbins,d0Cuts);
    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    if(!fIsMC)
      analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    else
      analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUsePID(true);
    if(fUseStrongPID) {
      analysiscuts->SetPidOption(1); //0=kConservative,1=kStrong
      analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
    }
    else analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb6080);

    analysiscuts->SetOptPileup(false);
    analysiscuts->SetMinCentrality(mincen);
    analysiscuts->SetMaxCentrality(maxcen);

    analysiscuts->SetMinPtCandidate(ptmin);
    analysiscuts->SetMaxPtCandidate(ptmax);
    analysiscuts->SetRemoveDaughtersFromPrim(false);

    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;

    analysiscuts->PrintAll();
    TString triggername = "kINT7";
    if(fIsMC) triggername = "kMB";
    TString pidname = "";
    if(fUseStrongPID) pidname = Form("_strongPIDpt%0.f", maxPtstrongPID);

    TFile* fout=new TFile(Form("DstoKKpiCuts_6080_loose%s_Raa_%s_pt%0.f_%0.f.root", pidname.Data(), triggername.Data(), ptmin, ptmax),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

AliRDHFCutsDstoKKpi* MakeFileForCutsDs6080_FiltTreeCreator2018(bool fUseStrongPID = false, double maxPtstrongPID = 8.0, bool fIsMC=false, double ptmin=2., double ptmax=50.) {

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3,1.e10);

    float mincen=60.;
    float maxcen=80;

    const int nptbins=2;
    float ptbins[nptbins+1]={0.,5.,50.};

    const int nvars = 20;

    float** anacutsval=new float*[nvars];
    for(int ic=0;ic<nvars;ic++){anacutsval[ic]=new float[nptbins];}

    anacutsval[0][0]=0.25;
    anacutsval[1][0]=0.4;
    anacutsval[2][0]=0.4;
    anacutsval[3][0]=0.;
    anacutsval[4][0]=0.;
    anacutsval[5][0]=0.;
    anacutsval[6][0]=0.05;
    anacutsval[7][0]=0.02;
    anacutsval[8][0]=0.;
    anacutsval[9][0]=0.95;
    anacutsval[10][0]=0.;
    anacutsval[11][0]=1000.0;
    anacutsval[12][0]=0.010;
    anacutsval[13][0]=0.001;
    anacutsval[14][0]=0.;
    anacutsval[15][0]=1.;
    anacutsval[16][0]=0.;
    anacutsval[17][0]=0.;
    anacutsval[18][0]=2.;
    anacutsval[19][0]=0.95;

    anacutsval[0][1]=0.3;
    anacutsval[1][1]=0.4;
    anacutsval[2][1]=0.4;
    anacutsval[3][1]=0.;
    anacutsval[4][1]=0.;
    anacutsval[5][1]=0.;
    anacutsval[6][1]=0.06;
    anacutsval[7][1]=0.02;
    anacutsval[8][1]=0.;
    anacutsval[9][1]=0.92;
    anacutsval[10][1]=0.;
    anacutsval[11][1]=100000.;
    anacutsval[12][1]=0.015;
    anacutsval[13][1]=0.0001;
    anacutsval[14][1]=-1.;
    anacutsval[15][1]=1.;
    anacutsval[16][1]=0.;
    anacutsval[17][1]=0.;
    anacutsval[18][1]=2.;
    anacutsval[19][1]=-1.;

    AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Ds Analysis and CF");
    analysiscuts->SetUsePreSelect(1);
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);

    TString cent="";

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    if(!fIsMC) analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    else analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUsePID(true);
    if(fUseStrongPID) {
      analysiscuts->SetPidOption(1); //0=kConservative,1=kStrong
      analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
    }
    else analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb6080);

    analysiscuts->SetOptPileup(false);
    if(!fIsMC) {
      analysiscuts->SetMinCentrality(mincen);
      analysiscuts->SetMaxCentrality(maxcen);
    }
    else {
      analysiscuts->SetMinCentrality(0);
      analysiscuts->SetMaxCentrality(100);
    }
    cent=Form("%.0f%.0f",mincen,maxcen);

    analysiscuts->SetMinPtCandidate(ptmin);
    analysiscuts->SetMaxPtCandidate(ptmax);
    analysiscuts->SetRemoveDaughtersFromPrim(false);

    cout<<"This is the object I'm going to save:"<<nptbins<<endl;

    analysiscuts->PrintAll();
    TString triggername = "kINT7";
    if(fIsMC) triggername = "kMB";
    TString pidname = "";
    if(fUseStrongPID) pidname = Form("_strongPIDpt%0.f", maxPtstrongPID);

    TFile* fout=new TFile(Form("DstoKKpiCuts_6080_filttreecreator%s_Raa_%s_pt%0.f_%0.f.root",pidname.Data(),triggername.Data(),ptmin,ptmax),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}
