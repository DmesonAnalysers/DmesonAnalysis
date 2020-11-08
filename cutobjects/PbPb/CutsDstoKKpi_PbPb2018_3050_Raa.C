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
// 1) MakeFileForCutsDs3050_Central2015 --> central cuts of 2015 analysis
// 2) MakeFileForCutsDs3050_Central2018 --> central cuts of 2018 analysis
// 3) MakeFileForCutsDs3050_Loose2018 --> loose cuts of 2018 analysis (sparse)
// 4) MakeFileForCutsDs3050_FiltTreeCreator2018 --> filtering cuts of 2018 analysis, locally on aliceml, SQM  (tree creator)
// 5) MakeFileForCutsDs3050_FiltTreeCreator2018QM --> filtering cuts of 2018 analysis, application on grid, QM (tree creator)
// 6) MakeFileForCutsDs3050_Central2018_Pass3 --> central cuts of 2018 analysis + Pass3 updates
// 7) MakeFileForCutsDs3050_Filt2018_Pass3 --> filtering cuts of 2018 analysis, application on grid + Pass3 updates
//____________________________________________________________________________________________________//

AliRDHFCutsDstoKKpi* MakeFileForCutsDs3050_Central2015(bool fUseStrongPID = true, double maxPtstrongPID = 8.0, bool fIsMC=false, double ptmin=2., double ptmax=50.) {

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3,1.e10);

    float mincen=30.;
    float maxcen=50.;

    const int nptbins=6;
    float* ptbins;
    ptbins=new float[nptbins+1];
    ptbins[0]=2.;
    ptbins[1]=4.;
    ptbins[2]=6.;
    ptbins[3]=8.;
    ptbins[4]=12.;
    ptbins[5]=16.;
    ptbins[6]=24.;

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
        anacutsval[6][ipt]=0.06;
        anacutsval[7][ipt]=0.02;
        anacutsval[8][ipt]=0.;
        anacutsval[9][ipt]=0.94;
        anacutsval[10][ipt]=0.;
        anacutsval[11][ipt]=1000.0;
        anacutsval[12][ipt]=0.015;
        anacutsval[13][ipt]=0.001;
        anacutsval[14][ipt]=0.;
        anacutsval[15][ipt]=1.;
        anacutsval[16][ipt]=0.;
        anacutsval[17][ipt]=0.;
        anacutsval[18][ipt]=0.;
        anacutsval[19][ipt]=0.94;

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

    anacutsval[6][0]=0.02;   //sigmavert
    anacutsval[7][0]=0.04;   //decay length
    anacutsval[9][0]=0.997;   //cosP
    anacutsval[12][0]=0.004; //Mass Phi
    anacutsval[14][0]=0.2;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][0]=0.8;  //CosP labFrame
    anacutsval[16][0]=0.03;  //decayXY
    anacutsval[17][0]=0.;    //normdecay
    anacutsval[18][0]=9.;   //normdecayXY
    anacutsval[19][0]=0.997;  //CosPXY

    anacutsval[6][1]=0.02;   //sigmavert
    anacutsval[7][1]=0.05;   //decay length
    anacutsval[9][1]=0.99;   //cosP
    anacutsval[12][1]=0.005; //Mass Phi
    anacutsval[14][1]=0.15;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][1]=0.7;  //CosP labFrame
    anacutsval[16][1]=0.05;  //decayXY
    anacutsval[17][1]=0.;    //normdecay
    anacutsval[18][1]=8.0;   //normdecayXY
    anacutsval[19][1]=0.99;  //CosPXY

    anacutsval[6][2]=0.025;   //sigmavert
    anacutsval[7][2]=0.05;   //decay length
    anacutsval[9][2]=0.98;   //cosP
    anacutsval[12][2]=0.006; //Mass Phi
    anacutsval[14][2]=0.05;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][2]=0.85;  //CosP labFrame
    anacutsval[16][2]=0.05;  //decayXY
    anacutsval[17][2]=0.;    //normdecay
    anacutsval[18][2]=7.0;   //normdecayXY
    anacutsval[19][2]=0.98;  //CosPXY

    anacutsval[6][3]=0.025;   //sigmavert
    anacutsval[7][3]=0.05;   //decay length
    anacutsval[9][3]=0.98;   //cosP
    anacutsval[12][3]=0.005; //Mass Phi
    anacutsval[14][3]=0.15;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][3]=0.85;  //CosP labFrame
    anacutsval[16][3]=0.05;  //decayXY
    anacutsval[17][3]=0.;    //normdecay
    anacutsval[18][3]=6.0;   //normdecayXY
    anacutsval[19][3]=0.98;  //CosPXY

    anacutsval[6][4]=0.02;   //sigmavert
    anacutsval[7][4]=0.04;   //decay length
    anacutsval[9][4]=0.98;   //cosP
    anacutsval[12][4]=0.005; //Mass Phi
    anacutsval[14][4]=0.0;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][4]=1.0;  //CosP labFrame
    anacutsval[16][4]=0.04;  //decayXY
    anacutsval[17][4]=0.;    //normdecay
    anacutsval[18][4]=5.0;   //normdecayXY
    anacutsval[19][4]=0.98;  //CosPXY

    anacutsval[6][5]=0.02;   //sigmavert
    anacutsval[7][5]=0.04;   //decay length
    anacutsval[9][5]=0.98;   //cosP
    anacutsval[12][5]=0.005; //Mass Phi
    anacutsval[14][5]=0.0;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][5]=1.0;  //CosP labFrame
    anacutsval[16][5]=0.04;  //decayXY
    anacutsval[17][5]=0.;    //normdecay
    anacutsval[18][5]=5.0;   //normdecayXY
    anacutsval[19][5]=0.98;  //CosPXY

    float topomCuts[nptbins] = {1.5,1.5,2.5,2.5,2.5,2.5}; //topomatic

    AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Ds Analysis and CF");
    analysiscuts->SetUsePreSelect(1);
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->Setd0MeasMinusExpCut(nptbins,topomCuts);
    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    if(!fIsMC)
      analysiscuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kSemiCentral);
    else
      analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUsePID(true);
    if(fUseStrongPID) {
      analysiscuts->SetPidOption(1); //0=kConservative,1=kStrong
      analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
    }
    else analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb3050);

    analysiscuts->SetOptPileup(false);
    analysiscuts->SetMinCentrality(mincen);
    analysiscuts->SetMaxCentrality(maxcen);

    analysiscuts->SetMinPtCandidate(ptmin);
    analysiscuts->SetMaxPtCandidate(ptmax);
    analysiscuts->SetRemoveDaughtersFromPrim(false);

    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;

    analysiscuts->PrintAll();
    TString triggername = "kINT7_kSemiCentral";
    if(fIsMC) triggername = "kMB";
    TString pidname = "";
    if(fUseStrongPID) pidname = Form("_strongPIDpt%0.f", maxPtstrongPID);

    TFile* fout=new TFile(Form("DstoKKpiCuts2015_3050_central%s_Raa_%s_pt%0.f_%0.f.root", pidname.Data(), triggername.Data(), ptmin, ptmax),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

AliRDHFCutsDstoKKpi* MakeFileForCutsDs3050_Central2018(bool fUseStrongPID = true, double maxPtstrongPID = 8.0, bool fIsMC = false, int addTrackCut = 0, double ptmin=2., double ptmax=50.) {

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetMinNCrossedRowsTPC(70);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3,1.e10);

    float mincen=30.;
    float maxcen=50.;

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

    anacutsval[6][0]=0.020;   //sigmavert
    anacutsval[7][0]=0.04;   //decay length
    anacutsval[9][0]=0.985;   //cosP
    anacutsval[12][0]=0.005; //Mass Phi
    anacutsval[14][0]=0.2;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][0]=1.0;  //CosP labFrame
    anacutsval[16][0]=0.04;  //decayXY
    anacutsval[17][0]=0.;    //normdecay
    anacutsval[18][0]=8.0;   //normdecayXY
    anacutsval[19][0]=0.99;  //CosPXY

    anacutsval[6][1]=0.025;   //sigmavert
    anacutsval[7][1]=0.04;   //decay length
    anacutsval[9][1]=0.985;   //cosP
    anacutsval[12][1]=0.005; //Mass Phi
    anacutsval[14][1]=0.15;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][1]=1.0;  //CosP labFrame
    anacutsval[16][1]=0.04;  //decayXY
    anacutsval[17][1]=0.;    //normdecay
    anacutsval[18][1]=8.0;   //normdecayXY
    anacutsval[19][1]=0.99;  //CosPXY

    anacutsval[6][2]=0.030;   //sigmavert
    anacutsval[7][2]=0.04;   //decay length
    anacutsval[9][2]=0.985;   //cosP
    anacutsval[12][2]=0.005; //Mass Phi
    anacutsval[14][2]=0.15;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][2]=1.0;   //CosP labFrame
    anacutsval[16][2]=0.04;  //decayXY
    anacutsval[17][2]=0.;    //normdecay
    anacutsval[18][2]=7.0;   //normdecayXY
    anacutsval[19][2]=0.99;  //CosPXY

    anacutsval[6][3]=0.030;   //sigmavert
    anacutsval[7][3]=0.04;   //decay length
    anacutsval[9][3]=0.98;   //cosP
    anacutsval[12][3]=0.005; //Mass Phi
    anacutsval[14][3]=0.15;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][3]=1.0;  //CosP labFrame
    anacutsval[16][3]=0.04;   //decayXY
    anacutsval[17][3]=0.;    //normdecay
    anacutsval[18][3]=7.0;   //normdecayXY
    anacutsval[19][3]=0.985;  //CosPXY

    anacutsval[6][4]=0.030;   //sigmavert
    anacutsval[7][4]=0.04;   //decay length
    anacutsval[9][4]=0.98;   //cosP
    anacutsval[12][4]=0.005; //Mass Phi
    anacutsval[14][4]=0.05;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][4]=1.0;  //CosP labFrame
    anacutsval[16][4]=0.04;   //decayXY
    anacutsval[17][4]=0.;    //normdecay
    anacutsval[18][4]=7.0;   //normdecayXY
    anacutsval[19][4]=0.985;  //CosPXY

    anacutsval[6][5]=0.030;   //sigmavert
    anacutsval[7][5]=0.05;   //decay length
    anacutsval[9][5]=0.98;   //cosP
    anacutsval[12][5]=0.005; //Mass Phi
    anacutsval[14][5]=0.05;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][5]=1.0;  //CosP labFrame
    anacutsval[16][5]=0.05;   //decayXY
    anacutsval[17][5]=0.;    //normdecay
    anacutsval[18][5]=6.0;   //normdecayXY
    anacutsval[19][5]=0.985;  //CosPXY

    anacutsval[6][6]=0.035;   //sigmavert
    anacutsval[7][6]=0.05;   //decay length
    anacutsval[9][6]=0.98;   //cosP
    anacutsval[12][6]=0.005; //Mass Phi
    anacutsval[14][6]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][6]=1.0;  //CosP labFrame
    anacutsval[16][6]=0.05;   //decayXY
    anacutsval[17][6]=0.;    //normdecay
    anacutsval[18][6]=6.0;   //normdecayXY
    anacutsval[19][6]=0.985;  //CosPXY

    anacutsval[6][7]=0.040;   //sigmavert
    anacutsval[7][7]=0.05;   //decay length
    anacutsval[9][7]=0.975;   //cosP
    anacutsval[12][7]=0.005; //Mass Phi
    anacutsval[14][7]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][7]=1.0;  //CosP labFrame
    anacutsval[16][7]=0.05;   //decayXY
    anacutsval[17][7]=0.;    //normdecay
    anacutsval[18][7]=5.0;   //normdecayXY
    anacutsval[19][7]=0.98;  //CosPXY

    anacutsval[6][8]=0.045;   //sigmavert
    anacutsval[7][8]=0.05;   //decay length
    anacutsval[9][8]=0.97;   //cosP
    anacutsval[12][8]=0.005; //Mass Phi
    anacutsval[14][8]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][8]=1.0;  //CosP labFrame
    anacutsval[16][8]=0.05;   //decayXY
    anacutsval[17][8]=0.;    //normdecay
    anacutsval[18][8]=4.0;   //normdecayXY
    anacutsval[19][8]=0.975;  //CosPXY

    float topomCuts[nptbins] = {1.5,2.,2.,2.5,2.5,2.5,3.,3.,3.}; //topomatic
    float d0Cuts[nptbins] = {70.,70.,70.,70.,70.,70.,70.,70.,70.}; //impparXY

    AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Ds Analysis and CF");
    analysiscuts->SetUsePreSelect(1);
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetMinNumTPCClsForPID(50);
    analysiscuts->Setd0MeasMinusExpCut(nptbins,topomCuts);
    analysiscuts->Setd0Cut(nptbins,d0Cuts);
    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    if(!fIsMC)
      analysiscuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kSemiCentral);
    else
      analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUsePID(true);
    if(fUseStrongPID) {
      analysiscuts->SetPidOption(1); //0=kConservative,1=kStrong
      analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
    }
    else analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb3050);

    analysiscuts->SetOptPileup(false);
    if(!fIsMC) {
      analysiscuts->SetMinCentrality(mincen);
      analysiscuts->SetMaxCentrality(maxcen);
    }
    else {
      analysiscuts->SetMinCentrality(0);
      analysiscuts->SetMaxCentrality(100);
    }

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
    TString triggername = "kINT7_kSemiCentral";
    if(fIsMC) triggername = "kMB";
    TString pidname = "";
    if(fUseStrongPID) pidname = Form("_strongPIDpt%0.f", maxPtstrongPID);

    TFile* fout=new TFile(Form("DstoKKpiCuts_3050_central%s_Raa_%s%s_pt%0.f_%0.f_pass1.root", pidname.Data(), triggername.Data(), trackCutName.Data(), ptmin, ptmax),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

AliRDHFCutsDstoKKpi* MakeFileForCutsDs3050_Loose2018(bool fUseStrongPID = true, double maxPtstrongPID = 8.0, bool fIsMC=false, double ptmin=2., double ptmax=50.) {

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3,1.e10);

    float mincen=30.;
    float maxcen=50.;

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
    anacutsval[7][0]=0.03;   //decay length
    anacutsval[9][0]=0.975;   //cosP
    anacutsval[12][0]=0.008; //Mass Phi
    anacutsval[14][0]=0.1;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][0]=1.0;  //CosP labFrame
    anacutsval[16][0]=0.;  //decayXY
    anacutsval[17][0]=0.;    //normdecay
    anacutsval[18][0]=5.0;   //normdecayXY
    anacutsval[19][0]=0.98;  //CosPXY

    anacutsval[6][1]=0.045;   //sigmavert
    anacutsval[7][1]=0.03;   //decay length
    anacutsval[9][1]=0.975;   //cosP
    anacutsval[12][1]=0.008; //Mass Phi
    anacutsval[14][1]=0.05;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][1]=1.0;  //CosP labFrame
    anacutsval[16][1]=0.;  //decayXY
    anacutsval[17][1]=0.;    //normdecay
    anacutsval[18][1]=5.0;   //normdecayXY
    anacutsval[19][1]=0.98;  //CosPXY

    anacutsval[6][2]=0.045;   //sigmavert
    anacutsval[7][2]=0.03;   //decay length
    anacutsval[9][2]=0.975;   //cosP
    anacutsval[12][2]=0.008; //Mass Phi
    anacutsval[14][2]=0.05;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][2]=1.0;   //CosP labFrame
    anacutsval[16][2]=0.;  //decayXY
    anacutsval[17][2]=0.;    //normdecay
    anacutsval[18][2]=5.0;   //normdecayXY
    anacutsval[19][2]=0.98;  //CosPXY

    anacutsval[6][3]=0.045;   //sigmavert
    anacutsval[7][3]=0.03;   //decay length
    anacutsval[9][3]=0.97;   //cosP
    anacutsval[12][3]=0.008; //Mass Phi
    anacutsval[14][3]=0.05;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][3]=1.0;  //CosP labFrame
    anacutsval[16][3]=0.;   //decayXY
    anacutsval[17][3]=0.;    //normdecay
    anacutsval[18][3]=5.0;   //normdecayXY
    anacutsval[19][3]=0.975;  //CosPXY

    anacutsval[6][4]=0.045;   //sigmavert
    anacutsval[7][4]=0.03;   //decay length
    anacutsval[9][4]=0.97;   //cosP
    anacutsval[12][4]=0.008; //Mass Phi
    anacutsval[14][4]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][4]=1.0;  //CosP labFrame
    anacutsval[16][4]=0.;   //decayXY
    anacutsval[17][4]=0.;    //normdecay
    anacutsval[18][4]=4.0;   //normdecayXY
    anacutsval[19][4]=0.975;  //CosPXY

    anacutsval[6][5]=0.045;   //sigmavert
    anacutsval[7][5]=0.04;   //decay length
    anacutsval[9][5]=0.965;   //cosP
    anacutsval[12][5]=0.008; //Mass Phi
    anacutsval[14][5]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][5]=1.0;  //CosP labFrame
    anacutsval[16][5]=0.;   //decayXY
    anacutsval[17][5]=0.;    //normdecay
    anacutsval[18][5]=4.0;   //normdecayXY
    anacutsval[19][5]=0.97;  //CosPXY

    anacutsval[6][6]=0.045;   //sigmavert
    anacutsval[7][6]=0.04;   //decay length
    anacutsval[9][6]=0.965;   //cosP
    anacutsval[12][6]=0.008; //Mass Phi
    anacutsval[14][6]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][6]=1.0;  //CosP labFrame
    anacutsval[16][6]=0.;   //decayXY
    anacutsval[17][6]=0.;    //normdecay
    anacutsval[18][6]=3.0;   //normdecayXY
    anacutsval[19][6]=0.97;  //CosPXY

    anacutsval[6][7]=0.055;   //sigmavert
    anacutsval[7][7]=0.04;   //decay length
    anacutsval[9][7]=0.96;   //cosP
    anacutsval[12][7]=0.008; //Mass Phi
    anacutsval[14][7]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][7]=1.0;  //CosP labFrame
    anacutsval[16][7]=0.;   //decayXY
    anacutsval[17][7]=0.;    //normdecay
    anacutsval[18][7]=3.0;   //normdecayXY
    anacutsval[19][7]=0.965;  //CosPXY

    anacutsval[6][8]=0.055;   //sigmavert
    anacutsval[7][8]=0.04;   //decay length
    anacutsval[9][8]=0.96;   //cosP
    anacutsval[12][8]=0.008; //Mass Phi
    anacutsval[14][8]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][8]=1.0;  //CosP labFrame
    anacutsval[16][8]=0.;   //decayXY
    anacutsval[17][8]=0.;    //normdecay
    anacutsval[18][8]=3.0;   //normdecayXY
    anacutsval[19][8]=0.965;  //CosPXY

    float topomCuts[nptbins] = {3.,4.,4.,5.,5.,5.,6.,6.,6.}; //topomatic
    float d0Cuts[nptbins] = {300.,300.,300.,300.,300.,300.,300.,300.,300.}; //impparXY

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
      analysiscuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kSemiCentral);
    else
      analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUsePID(true);
    if(fUseStrongPID) {
      analysiscuts->SetPidOption(1); //0=kConservative,1=kStrong
      analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
    }
    else analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb3050);

    analysiscuts->SetOptPileup(false);
    if(!fIsMC) {
      analysiscuts->SetMinCentrality(mincen);
      analysiscuts->SetMaxCentrality(maxcen);
    }
    else {
      analysiscuts->SetMinCentrality(0);
      analysiscuts->SetMaxCentrality(100);
    }

    analysiscuts->SetMinPtCandidate(ptmin);
    analysiscuts->SetMaxPtCandidate(ptmax);
    analysiscuts->SetRemoveDaughtersFromPrim(false);

    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;

    analysiscuts->PrintAll();
    TString triggername = "kINT7_kSemiCentral";
    if(fIsMC) triggername = "kMB";
    TString pidname = "";
    if(fUseStrongPID) pidname = Form("_strongPIDpt%0.f", maxPtstrongPID);

    TFile* fout=new TFile(Form("DstoKKpiCuts2018_3050_loose%s_Raa_%s_pt%0.f_%0.f.root", pidname.Data(), triggername.Data(), ptmin, ptmax),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

AliRDHFCutsDstoKKpi* MakeFileForCutsDs3050_FiltTreeCreator2018(bool fUseStrongPID = false, double maxPtstrongPID = 8.0, bool fIsMC=false, double ptmin=2., double ptmax=50.) {

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    //default
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    esdTrackCuts->SetEtaRange(-0.8,0.8);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.4,1.e10);
    esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0025*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");

    float mincen=30.;
    float maxcen=50;

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
    anacutsval[6][0]=0.04;
    anacutsval[7][0]=0.03;
    anacutsval[8][0]=0.;
    anacutsval[9][0]=0.96;
    anacutsval[10][0]=0.;
    anacutsval[11][0]=1000.0;
    anacutsval[12][0]=0.010;
    anacutsval[13][0]=0.001;
    anacutsval[14][0]=0.;
    anacutsval[15][0]=1.;
    anacutsval[16][0]=0.;
    anacutsval[17][0]=0.;
    anacutsval[18][0]=3.;
    anacutsval[19][0]=0.96;

    anacutsval[0][1]=0.3;
    anacutsval[1][1]=0.4;
    anacutsval[2][1]=0.4;
    anacutsval[3][1]=0.;
    anacutsval[4][1]=0.;
    anacutsval[5][1]=0.;
    anacutsval[6][1]=0.06;
    anacutsval[7][1]=0.03;
    anacutsval[8][1]=0.;
    anacutsval[9][1]=0.95;
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

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    if(!fIsMC) analysiscuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kSemiCentral);
    else analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUsePID(true);
    if(fUseStrongPID) {
      analysiscuts->SetPidOption(1); //0=kConservative,1=kStrong
      analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
    }
    else analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb3050);

    analysiscuts->SetOptPileup(false);
    if(!fIsMC) {
      analysiscuts->SetMinCentrality(mincen);
      analysiscuts->SetMaxCentrality(maxcen);
    }
    else {
      analysiscuts->SetMinCentrality(0);
      analysiscuts->SetMaxCentrality(100);
    }

    analysiscuts->SetMinPtCandidate(ptmin);
    analysiscuts->SetMaxPtCandidate(ptmax);
    analysiscuts->SetRemoveDaughtersFromPrim(false);

    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;

    analysiscuts->PrintAll();
    TString triggername = "kINT7_kSemiCentral";
    if(fIsMC) triggername = "kMB";
    TString pidname = "";
    if(fUseStrongPID) pidname = Form("_strongPIDpt%0.f", maxPtstrongPID);

    TFile* fout=new TFile(Form("DstoKKpiCuts_3050_filttreecreator%s_Raa_%s_pt%0.f_%0.f.root",pidname.Data(),triggername.Data(),ptmin,ptmax),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

AliRDHFCutsDstoKKpi* MakeFileForCutsDs3050_FiltTreeCreator2018QM(bool fIsMC=false, double ptmin=2., double ptmax=50.) {

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.4,1.e10);

    float mincen=30.;
    float maxcen=50;

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
    anacutsval[6][0]=0.045;
    anacutsval[7][0]=0.02;
    anacutsval[8][0]=0.;
    anacutsval[9][0]=0.96;
    anacutsval[10][0]=0.;
    anacutsval[11][0]=1000.0;
    anacutsval[12][0]=0.015;
    anacutsval[13][0]=0.001;
    anacutsval[14][0]=0.;
    anacutsval[15][0]=1.;
    anacutsval[16][0]=0.;
    anacutsval[17][0]=0.;
    anacutsval[18][0]=2.;
    anacutsval[19][0]=0.96;

    anacutsval[0][1]=0.3;
    anacutsval[1][1]=0.4;
    anacutsval[2][1]=0.4;
    anacutsval[3][1]=0.;
    anacutsval[4][1]=0.;
    anacutsval[5][1]=0.;
    anacutsval[6][1]=0.06;
    anacutsval[7][1]=0.02;
    anacutsval[8][1]=0.;
    anacutsval[9][1]=0.9;
    anacutsval[10][1]=0.;
    anacutsval[11][1]=1000.;
    anacutsval[12][1]=0.015;
    anacutsval[13][1]=0.0001;
    anacutsval[14][1]=0.;
    anacutsval[15][1]=1.;
    anacutsval[16][1]=0.;
    anacutsval[17][1]=0.;
    anacutsval[18][1]=2.;
    anacutsval[19][1]=0.9;

    AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Ds Analysis and CF");
    analysiscuts->SetUsePreSelect(1);
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    if(!fIsMC)
      analysiscuts->SetUseTimeRangeCutForPbPb2018(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    if(!fIsMC) analysiscuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kSemiCentral);
    else analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUsePID(true);
    analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb3050);

    analysiscuts->SetOptPileup(false);
    if(!fIsMC) {
      analysiscuts->SetMinCentrality(mincen);
      analysiscuts->SetMaxCentrality(maxcen);
    }
    else {
      analysiscuts->SetMinCentrality(0);
      analysiscuts->SetMaxCentrality(100);
    }

    analysiscuts->SetMinPtCandidate(ptmin);
    analysiscuts->SetMaxPtCandidate(ptmax);
    analysiscuts->SetRemoveDaughtersFromPrim(false);

    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;

    analysiscuts->PrintAll();
    TString triggername = "kINT7_kSemiCentral";
    if(fIsMC) triggername = "kMB";

    TFile* fout=new TFile(Form("DstoKKpiCuts_3050_filttreecreatorQM_Raa_%s_pt%0.f_%0.f.root",triggername.Data(),ptmin,ptmax),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

AliRDHFCutsDstoKKpi* MakeFileForCutsDs3050_Central2018_Pass3(bool fUseStrongPID = true, double maxPtstrongPID = 8.0, bool fIsMC = false, int addTrackCut = 0, double ptmin=2., double ptmax=50.) {

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMaxChi2PerClusterTPC(2.5);
    esdTrackCuts->SetMinNCrossedRowsTPC(70);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.4,1.e10);

    float mincen=30.;
    float maxcen=50.;

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
        anacutsval[1][ipt]=0.4;
        anacutsval[2][ipt]=0.4;
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

    anacutsval[6][0]=0.020;   //sigmavert
    anacutsval[7][0]=0.04;   //decay length
    anacutsval[9][0]=0.985;   //cosP
    anacutsval[12][0]=0.005; //Mass Phi
    anacutsval[14][0]=0.2;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][0]=1.0;  //CosP labFrame
    anacutsval[16][0]=0.04;  //decayXY
    anacutsval[17][0]=0.;    //normdecay
    anacutsval[18][0]=8.0;   //normdecayXY
    anacutsval[19][0]=0.99;  //CosPXY

    anacutsval[6][1]=0.025;   //sigmavert
    anacutsval[7][1]=0.04;   //decay length
    anacutsval[9][1]=0.985;   //cosP
    anacutsval[12][1]=0.005; //Mass Phi
    anacutsval[14][1]=0.15;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][1]=1.0;  //CosP labFrame
    anacutsval[16][1]=0.04;  //decayXY
    anacutsval[17][1]=0.;    //normdecay
    anacutsval[18][1]=8.0;   //normdecayXY
    anacutsval[19][1]=0.99;  //CosPXY

    anacutsval[6][2]=0.030;   //sigmavert
    anacutsval[7][2]=0.04;   //decay length
    anacutsval[9][2]=0.985;   //cosP
    anacutsval[12][2]=0.005; //Mass Phi
    anacutsval[14][2]=0.15;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][2]=1.0;   //CosP labFrame
    anacutsval[16][2]=0.04;  //decayXY
    anacutsval[17][2]=0.;    //normdecay
    anacutsval[18][2]=7.0;   //normdecayXY
    anacutsval[19][2]=0.99;  //CosPXY

    anacutsval[6][3]=0.030;   //sigmavert
    anacutsval[7][3]=0.04;   //decay length
    anacutsval[9][3]=0.98;   //cosP
    anacutsval[12][3]=0.005; //Mass Phi
    anacutsval[14][3]=0.15;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][3]=1.0;  //CosP labFrame
    anacutsval[16][3]=0.04;   //decayXY
    anacutsval[17][3]=0.;    //normdecay
    anacutsval[18][3]=7.0;   //normdecayXY
    anacutsval[19][3]=0.985;  //CosPXY

    anacutsval[6][4]=0.030;   //sigmavert
    anacutsval[7][4]=0.04;   //decay length
    anacutsval[9][4]=0.98;   //cosP
    anacutsval[12][4]=0.005; //Mass Phi
    anacutsval[14][4]=0.05;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][4]=1.0;  //CosP labFrame
    anacutsval[16][4]=0.04;   //decayXY
    anacutsval[17][4]=0.;    //normdecay
    anacutsval[18][4]=7.0;   //normdecayXY
    anacutsval[19][4]=0.985;  //CosPXY

    anacutsval[6][5]=0.030;   //sigmavert
    anacutsval[7][5]=0.05;   //decay length
    anacutsval[9][5]=0.98;   //cosP
    anacutsval[12][5]=0.005; //Mass Phi
    anacutsval[14][5]=0.05;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][5]=1.0;  //CosP labFrame
    anacutsval[16][5]=0.05;   //decayXY
    anacutsval[17][5]=0.;    //normdecay
    anacutsval[18][5]=6.0;   //normdecayXY
    anacutsval[19][5]=0.985;  //CosPXY

    anacutsval[6][6]=0.035;   //sigmavert
    anacutsval[7][6]=0.05;   //decay length
    anacutsval[9][6]=0.98;   //cosP
    anacutsval[12][6]=0.005; //Mass Phi
    anacutsval[14][6]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][6]=1.0;  //CosP labFrame
    anacutsval[16][6]=0.05;   //decayXY
    anacutsval[17][6]=0.;    //normdecay
    anacutsval[18][6]=6.0;   //normdecayXY
    anacutsval[19][6]=0.985;  //CosPXY

    anacutsval[6][7]=0.040;   //sigmavert
    anacutsval[7][7]=0.05;   //decay length
    anacutsval[9][7]=0.975;   //cosP
    anacutsval[12][7]=0.005; //Mass Phi
    anacutsval[14][7]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][7]=1.0;  //CosP labFrame
    anacutsval[16][7]=0.05;   //decayXY
    anacutsval[17][7]=0.;    //normdecay
    anacutsval[18][7]=5.0;   //normdecayXY
    anacutsval[19][7]=0.98;  //CosPXY

    anacutsval[6][8]=0.045;   //sigmavert
    anacutsval[7][8]=0.05;   //decay length
    anacutsval[9][8]=0.97;   //cosP
    anacutsval[12][8]=0.005; //Mass Phi
    anacutsval[14][8]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][8]=1.0;  //CosP labFrame
    anacutsval[16][8]=0.05;   //decayXY
    anacutsval[17][8]=0.;    //normdecay
    anacutsval[18][8]=4.0;   //normdecayXY
    anacutsval[19][8]=0.975;  //CosPXY

    float topomCuts[nptbins] = {1.5,2.,2.,2.5,2.5,2.5,3.,3.,3.}; //topomatic
    float d0Cuts[nptbins] = {70.,70.,70.,70.,70.,70.,70.,70.,70.}; //impparXY

    AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Ds Analysis and CF");
    analysiscuts->SetUsePreSelect(1);
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetMinNumTPCClsForPID(50);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->Setd0MeasMinusExpCut(nptbins,topomCuts);
    analysiscuts->Setd0Cut(nptbins,d0Cuts);
    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    if(!fIsMC)
      analysiscuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kSemiCentral);
    else
      analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUsePID(true);
    if(fUseStrongPID) {
      analysiscuts->SetPidOption(1); //0=kConservative,1=kStrong
      analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
    }
    else
      analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb3050);
    analysiscuts->SetOptPileup(false);
    if(!fIsMC) {
      analysiscuts->SetMinCentrality(mincen);
      analysiscuts->SetMaxCentrality(maxcen);
    }
    else {
      analysiscuts->SetMinCentrality(0);
      analysiscuts->SetMaxCentrality(100);
    }

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
    TString triggername = "kINT7_kSemiCentral";
    if(fIsMC) triggername = "kMB";
    TString pidname = "";
    if(fUseStrongPID) pidname = Form("_strongPIDpt%0.f", maxPtstrongPID);

    TFile* fout=new TFile(Form("DstoKKpiCuts_3050_central%s_Raa_%s%s_pt%0.f_%0.f_2018pass3.root", pidname.Data(), triggername.Data(), trackCutName.Data(), ptmin, ptmax),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

AliRDHFCutsDstoKKpi* MakeFileForCutsDs3050_Filt2018_Pass3(bool fIsMC=false, double ptmin=2., double ptmax=50.) {

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMaxChi2PerClusterTPC(2.5);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetMinNCrossedRowsTPC(70);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.4,1.e10);

    float mincen=30.;
    float maxcen=50;

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
    anacutsval[6][0]=0.045;
    anacutsval[7][0]=0.02;
    anacutsval[8][0]=0.;
    anacutsval[9][0]=0.96;
    anacutsval[10][0]=0.;
    anacutsval[11][0]=1000.0;
    anacutsval[12][0]=0.015;
    anacutsval[13][0]=0.001;
    anacutsval[14][0]=0.;
    anacutsval[15][0]=1.;
    anacutsval[16][0]=0.;
    anacutsval[17][0]=0.;
    anacutsval[18][0]=2.;
    anacutsval[19][0]=0.96;

    anacutsval[0][1]=0.3;
    anacutsval[1][1]=0.4;
    anacutsval[2][1]=0.4;
    anacutsval[3][1]=0.;
    anacutsval[4][1]=0.;
    anacutsval[5][1]=0.;
    anacutsval[6][1]=0.06;
    anacutsval[7][1]=0.02;
    anacutsval[8][1]=0.;
    anacutsval[9][1]=0.9;
    anacutsval[10][1]=0.;
    anacutsval[11][1]=1000.;
    anacutsval[12][1]=0.015;
    anacutsval[13][1]=0.0001;
    anacutsval[14][1]=0.;
    anacutsval[15][1]=1.;
    anacutsval[16][1]=0.;
    anacutsval[17][1]=0.;
    anacutsval[18][1]=2.;
    anacutsval[19][1]=0.9;

    AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Ds Analysis and CF");
    analysiscuts->SetUsePreSelect(1);
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetMinNumTPCClsForPID(50);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    if(!fIsMC)
      analysiscuts->SetUseTimeRangeCutForPbPb2018(true);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    if(!fIsMC) 
      analysiscuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kSemiCentral);
    else 
      analysiscuts->SetTriggerMask(AliVEvent::kMB);

    analysiscuts->SetUsePID(true);
    analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb3050);
    analysiscuts->SetOptPileup(false);
    if(!fIsMC) {
      analysiscuts->SetMinCentrality(mincen);
      analysiscuts->SetMaxCentrality(maxcen);
    }
    else {
      analysiscuts->SetMinCentrality(0);
      analysiscuts->SetMaxCentrality(100);
    }

    analysiscuts->SetMinPtCandidate(ptmin);
    analysiscuts->SetMaxPtCandidate(ptmax);
    analysiscuts->SetRemoveDaughtersFromPrim(false);

    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;

    analysiscuts->PrintAll();
    TString triggername = "kINT7_kSemiCentral";
    if(fIsMC) triggername = "kMB";

    TFile* fout=new TFile(Form("DstoKKpiCuts_3050_filt_Raa_%s_pt%0.f_%0.f_2018pass3.root",triggername.Data(),ptmin,ptmax),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}
