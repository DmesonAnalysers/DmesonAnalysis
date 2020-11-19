#if !defined (__CINT__) || defined (__CLING__)

#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>

#include "AliESDtrackCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"

#endif

//____________________________________________________________________________________________________//
// Methods:
// 1) MakeFileForCutsDplus010_Central2018 --> central cuts of 2015 and 2018 analyses
// 2) MakeFileForCutsDplus010_Loose2018 --> loose cuts of 2018 analysis (sparse)
// 3) MakeFileForCutsDplus010_FiltTreeCreator2018 --> filtering cuts of 2018 analysis (tree creator)
//____________________________________________________________________________________________________//

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi* MakeFileForCutsDplus010_Central2018(bool fUseStrongPID=true, double maxPtstrongPID=3.0, bool fIsMC=false, int addTrackCut = 0, int addGeoCut=0){

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetMaxChi2PerClusterTPC(2.5);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0060*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
    esdTrackCuts->SetPtRange(0.6,1.e10);

    float minc=0.,maxc=10.;
    const int nptbins=15;
    float* ptbins;
    ptbins=new float[nptbins+1];

    ptbins[0]=2.;
    ptbins[1]=3.;
    ptbins[2]=4.;
    ptbins[3]=5.;
    ptbins[4]=6.;
    ptbins[5]=7.;
    ptbins[6]=8.;
    ptbins[7]=9.;
    ptbins[8]=10.;
    ptbins[9]=11.;
    ptbins[10]=12.;
    ptbins[11]=14.;
    ptbins[12]=16.;
    ptbins[13]=24.;
    ptbins[14]=36.;
    ptbins[15]=50.;

    const int nvars=14;
    float** anacutsval;
    anacutsval=new float*[nvars];
    for(int ic=0;ic<nvars;ic++){anacutsval[ic]=new float[nptbins];}

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

    for(int ipt=0;ipt<nptbins;ipt++){
        anacutsval[0][ipt]=0.2; //minv
        anacutsval[1][ipt]=0.6; //ptK
        anacutsval[2][ipt]=0.6; //ptPi
        anacutsval[3][ipt]=0.; //d0K
        anacutsval[4][ipt]=0.; //d0Pi
        anacutsval[5][ipt]=0.; //dist12
        anacutsval[8][ipt]=0.0; //pM
        anacutsval[10][ipt]=0.0; //sumd02
        anacutsval[11][ipt]=10000000000.; //dca
    }

    int ic=6;//sigvert
    anacutsval[ic][0]=0.020;//2.0-3.0
    anacutsval[ic][1]=0.020;//3.0-4.0
    anacutsval[ic][2]=0.020;//4.0-5.0
    anacutsval[ic][3]=0.020;//5.0-6.0
    anacutsval[ic][4]=0.020;//6.0-7.0
    anacutsval[ic][5]=0.024;//7.0-8.0
    anacutsval[ic][6]=0.024;//8.0-9.0
    anacutsval[ic][7]=0.024;//9.0-10.0
    anacutsval[ic][8]=0.024;//10.0-11.0
    anacutsval[ic][9]=0.024;//11.0-12.0
    anacutsval[ic][10]=0.024;//12.0-14.0
    anacutsval[ic][11]=0.024;//14.0-16.0
    anacutsval[ic][12]=0.024;//16.0-24.0
    anacutsval[ic][13]=0.034;//24.0-36.0
    anacutsval[ic][14]=0.040;//36.0-50.0

    ic=7;//declen
    for(int ipt=0;ipt<nptbins;ipt++){
        anacutsval[ic][ipt]=0.14;
    }
    anacutsval[ic][0]=0.13;//2.0-3.0
    anacutsval[ic][1]=0.13;//3.0-4.0
    anacutsval[ic][2]=0.13;//4.0-5.0
    anacutsval[ic][3]=0.13;//5.0-6.0
    anacutsval[ic][10]=0.16;//12.0-14.0
    anacutsval[ic][11]=0.16;//14.0-16.0
    anacutsval[ic][12]=0.18;//16.0-24.0
    anacutsval[ic][13]=0.20;//24.0-36.0
    anacutsval[ic][14]=0.20;//36.0-50.0

    //cosp
    ic=9;
    anacutsval[ic][0]=0.998;//2.0-3.0
    anacutsval[ic][1]=0.998;//3.0-4.0
    anacutsval[ic][2]=0.998;//4.0-5.0
    anacutsval[ic][3]=0.997;//5.0-6.0
    anacutsval[ic][4]=0.997;//6.0-7.0
    anacutsval[ic][5]=0.995;//7.0-8.0
    anacutsval[ic][6]=0.992;//8.0-9.0
    anacutsval[ic][7]=0.992;//9.0-10.0
    anacutsval[ic][8]=0.992;//10.0-11.0
    anacutsval[ic][9]=0.992;//11.0-12.0
    anacutsval[ic][10]=0.990;//12.0-14.0
    anacutsval[ic][11]=0.990;//14.0-16.0
    anacutsval[ic][12]=0.990;//16.0-24.0
    anacutsval[ic][13]=0.990;//24.0-36.0
    anacutsval[ic][14]=0.990;//36.0-50.0

    ic=12;//ndlXY
    anacutsval[ic][0]=15.;//2.0-3.0
    anacutsval[ic][1]=14.;//3.0-4.0
    anacutsval[ic][2]=12.;//4.0-5.0
    anacutsval[ic][3]=12.;//5.0-6.0
    anacutsval[ic][4]=12.;//6.0-7.0
    anacutsval[ic][5]=12.;//7.0-8.0
    anacutsval[ic][6]=12.;//8.0-9.0
    anacutsval[ic][7]=12.;//9.0-10.0
    anacutsval[ic][8]=10.;//10.0-11.0
    anacutsval[ic][9]=10.;//11.0-12.0
    anacutsval[ic][10]=10.;//12.0-14.0
    anacutsval[ic][11]=10.;//14.0-16.0
    anacutsval[ic][12]=10.;//16.0-24.0
    anacutsval[ic][13]=10.;//24.0-36.0
    anacutsval[ic][14]=10.;//36.0-50.0

    ic=13;//cospXY
    anacutsval[ic][0]=0.999;//2.0-3.0
    anacutsval[ic][1]=0.999;//3.0-4.0
    anacutsval[ic][2]=0.999;//4.0-5.0
    anacutsval[ic][3]=0.998;//5.0-6.0
    anacutsval[ic][4]=0.998;//6.0-7.0
    anacutsval[ic][5]=0.996;//7.0-8.0
    anacutsval[ic][6]=0.994;//8.0-9.0
    anacutsval[ic][7]=0.994;//9.0-10.0
    anacutsval[ic][8]=0.994;//10.0-11.0
    anacutsval[ic][9]=0.994;//11.0-12.0
    anacutsval[ic][10]=0.990;//12.0-14.0
    anacutsval[ic][11]=0.990;//14.0-16.0
    anacutsval[ic][12]=0.990;//16.0-24.0
    anacutsval[ic][13]=0.990;//24.0-36.0
    anacutsval[ic][14]=0.990;//36.0-50.0

    float *d0cutsval=new float[nptbins];
    for(int ipt=0;ipt<nptbins;ipt++){ //d0
        d0cutsval[ipt]=60;
    }
    d0cutsval[0]=90;
    d0cutsval[12]=80;
    d0cutsval[13]=80;
    d0cutsval[14]=80;

    float *d0d0expcutsval=new float[nptbins];
    for(int ipt=0;ipt<nptbins;ipt++){ //d0d0exp
        d0d0expcutsval[ipt]=2.0;
    }
    d0d0expcutsval[13]=2.5;
    d0d0expcutsval[14]=2.5;

    AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->Setd0Cut(nptbins,d0cutsval);
    analysiscuts->Setd0MeasMinusExpCut(nptbins,d0d0expcutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetScaleNormDLxyBypOverPt(false);
    analysiscuts->SetMinNumTPCClsForPID(50);
    analysiscuts->SetUsePreSelect(1);

    analysiscuts->SetUsePID(true);
    AliAODPidHF* PidHF = NULL;
    if(fUseStrongPID) {
        analysiscuts->SetUseStrongPid(3);
        analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
        analysiscuts->SetMaxPStrongPidK(1);
        analysiscuts->SetMaxPStrongPidpi(1);
    }

    analysiscuts->SetUseImpParProdCorrCut(false);
    analysiscuts->SetOptPileup(false);
    if(!fIsMC) {
        analysiscuts->SetMinCentrality(minc);
        analysiscuts->SetMaxCentrality(maxc);
    }

    analysiscuts->SetRemoveTrackletOutliers(true);//added on June 28
    analysiscuts->SetCutOnzVertexSPD(3);//needed for Pb-Pb 2015

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kCentral);
    if(fIsMC) analysiscuts->SetTriggerMask(AliVEvent::kMB);

    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb010);

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

    TString geoCutName = "";
    switch(addGeoCut) {
        case 0:
            break;
        case 1:
            analysiscuts->SetUseCutGeoNcrNcl(true);
            geoCutName = "_geocutstd";
            break;
        case 2:
            analysiscuts->SetUseCutGeoNcrNcl(true);
            analysiscuts->ConfigureCutGeoNcrNcl(2, 130, 1.5, 0.85, 0.7);
            geoCutName = "_geocutDZ2Pt1dot5";
            break;
        case 3:
            analysiscuts->SetUseCutGeoNcrNcl(true);
            analysiscuts->ConfigureCutGeoNcrNcl(4, 130, 1.5, 0.85, 0.7);
            geoCutName = "_geocutDZ4Pt1dot5";
            break;
        case 4:
            analysiscuts->SetUseCutGeoNcrNcl(true);
            analysiscuts->ConfigureCutGeoNcrNcl(3, 130, 2., 0.85, 0.7);
            geoCutName = "_geocutDZ3Pt2";
            break;
        case 5:
            analysiscuts->SetUseCutGeoNcrNcl(true);
            analysiscuts->ConfigureCutGeoNcrNcl(3, 130, 4., 0.85, 0.7);
            geoCutName = "_geocutDZ3Pt4";
            break;
        case 6:
            analysiscuts->SetUseCutGeoNcrNcl(true);
            analysiscuts->ConfigureCutGeoNcrNcl(2, 130, 2., 0.85, 0.7);
            geoCutName = "_geocutDZ2Pt4";
            break;
    }

    analysiscuts->SetMinPtCandidate(2.);
    analysiscuts->SetMaxPtCandidate(50.);

    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;

    TString triggername="kINT7_kCentral";
    if(fIsMC) triggername="kMB";
    TString PIDsuffix="";
    if(fUseStrongPID) PIDsuffix=Form("_strongPIDpt%0.f",maxPtstrongPID);

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename=Form("DplustoKpipiCuts_010_central%s_Raa_%s%s%s.root",PIDsuffix.Data(),triggername.Data(),trackCutName.Data(),geoCutName.Data());
    TFile* fout=new TFile(filename.Data(),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi* MakeFileForCutsDplus010_Loose2018(bool fUseStrongPID=true, double maxPtstrongPID=3.0, bool fIsMC=false){

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMaxChi2PerClusterTPC(2.5);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0060*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
    esdTrackCuts->SetPtRange(0.6,1.e10);

    float minc=0.,maxc=10.;
    const int nptbins=15;
    float* ptbins;
    ptbins=new float[nptbins+1];

    ptbins[0]=2.;
    ptbins[1]=3.;
    ptbins[2]=4.;
    ptbins[3]=5.;
    ptbins[4]=6.;
    ptbins[5]=7.;
    ptbins[6]=8.;
    ptbins[7]=9.;
    ptbins[8]=10.;
    ptbins[9]=11.;
    ptbins[10]=12.;
    ptbins[11]=14.;
    ptbins[12]=16.;
    ptbins[13]=24.;
    ptbins[14]=36.;
    ptbins[15]=50.;

    const int nvars=14;
    float** anacutsval;
    anacutsval=new float*[nvars];
    for(int ic=0;ic<nvars;ic++){anacutsval[ic]=new float[nptbins];}

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

    for(int ipt=0;ipt<nptbins;ipt++){
        anacutsval[0][ipt]=0.2; //minv
        anacutsval[1][ipt]=0.6; //ptK
        anacutsval[2][ipt]=0.6; //ptPi
        anacutsval[3][ipt]=0.; //d0K
        anacutsval[4][ipt]=0.; //d0Pi
        anacutsval[5][ipt]=0.; //dist12
        anacutsval[8][ipt]=0.0; //pM
        anacutsval[10][ipt]=0.0; //sumd02
        anacutsval[11][ipt]=10000000000.; //dca
    }

    int ic=6;//sigvert
    anacutsval[ic][0]=0.024;//2.0-3.0
    anacutsval[ic][1]=0.024;//3.0-4.0
    anacutsval[ic][2]=0.024;//4.0-5.0
    anacutsval[ic][3]=0.026;//5.0-6.0
    anacutsval[ic][4]=0.026;//6.0-7.0
    anacutsval[ic][5]=0.026;//7.0-8.0
    anacutsval[ic][6]=0.026;//8.0-9.0
    anacutsval[ic][7]=0.026;//9.0-10.0
    anacutsval[ic][8]=0.026;//10.0-11.0
    anacutsval[ic][9]=0.026;//11.0-12.0
    anacutsval[ic][10]=0.032;//12.0-14.0
    anacutsval[ic][11]=0.032;//14.0-16.0
    anacutsval[ic][12]=0.036;//16.0-24.0
    anacutsval[ic][13]=0.050;//24.0-36.0
    anacutsval[ic][14]=0.050;//36.0-50.0

    ic=7;//declen
    for(int ipt=0;ipt<nptbins;ipt++){
        anacutsval[ic][ipt]=0.08;
    }
    anacutsval[ic][0]=0.06;//2.0-3.0
    anacutsval[ic][1]=0.08;//3.0-4.0
    anacutsval[ic][12]=0.10;//16.0-24.0
    anacutsval[ic][13]=0.10;//24.0-36.0
    anacutsval[ic][14]=0.10;//36.0-50.0

    //cosp
    ic=9;
    anacutsval[ic][0]=0.995;//2.0-3.0
    anacutsval[ic][1]=0.995;//3.0-4.0
    anacutsval[ic][2]=0.992;//4.0-5.0
    anacutsval[ic][3]=0.992;//5.0-6.0
    anacutsval[ic][4]=0.990;//6.0-7.0
    anacutsval[ic][5]=0.990;//7.0-8.0
    anacutsval[ic][6]=0.985;//8.0-9.0
    anacutsval[ic][7]=0.985;//9.0-10.0
    anacutsval[ic][8]=0.980;//10.0-11.0
    anacutsval[ic][9]=0.980;//11.0-12.0
    anacutsval[ic][10]=0.980;//12.0-14.0
    anacutsval[ic][11]=0.980;//14.0-16.0
    anacutsval[ic][12]=0.950;//16.0-24.0
    anacutsval[ic][13]=0.950;//24.0-36.0
    anacutsval[ic][14]=0.950;//36.0-50.0

    ic=12;//ndlXY
    anacutsval[ic][0]=10.;//2.0-3.0
    anacutsval[ic][1]=10.;//3.0-4.0
    anacutsval[ic][2]=10.;//4.0-5.0
    anacutsval[ic][3]=10.;//5.0-6.0
    anacutsval[ic][4]=8.;//6.0-7.0
    anacutsval[ic][5]=8.;//7.0-8.0
    anacutsval[ic][6]=8.;//8.0-9.0
    anacutsval[ic][7]=8.;//9.0-10.0
    anacutsval[ic][8]=8.;//10.0-11.0
    anacutsval[ic][9]=8.;//11.0-12.0
    anacutsval[ic][10]=6.;//12.0-14.0
    anacutsval[ic][11]=6.;//14.0-16.0
    anacutsval[ic][12]=6.;//16.0-24.0
    anacutsval[ic][13]=6.;//24.0-36.0
    anacutsval[ic][14]=6.;//36.0-50.0

    ic=13;//cospXY
    anacutsval[ic][0]=0.994;//2.0-3.0
    anacutsval[ic][1]=0.994;//3.0-4.0
    anacutsval[ic][2]=0.992;//4.0-5.0
    anacutsval[ic][3]=0.992;//5.0-6.0
    anacutsval[ic][4]=0.990;//6.0-7.0
    anacutsval[ic][5]=0.990;//7.0-8.0
    anacutsval[ic][6]=0.985;//8.0-9.0
    anacutsval[ic][7]=0.985;//9.0-10.0
    anacutsval[ic][8]=0.980;//10.0-11.0
    anacutsval[ic][9]=0.980;//11.0-12.0
    anacutsval[ic][10]=0.980;//12.0-14.0
    anacutsval[ic][11]=0.980;//14.0-16.0
    anacutsval[ic][12]=0.950;//16.0-24.0
    anacutsval[ic][13]=0.950;//24.0-36.0
    anacutsval[ic][14]=0.950;//36.0-50.0

    float *d0cutsval=new float[nptbins];
    for(int ipt=0;ipt<nptbins;ipt++){ //d0
        d0cutsval[ipt]=150;
    }

    float *d0d0expcutsval=new float[nptbins];
    for(int ipt=0;ipt<nptbins;ipt++){ //d0d0exp
        d0d0expcutsval[ipt]=4.5;
    }

    AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->Setd0Cut(nptbins,d0cutsval);
    analysiscuts->Setd0MeasMinusExpCut(nptbins,d0d0expcutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetMinNumTPCClsForPID(50);
    analysiscuts->SetScaleNormDLxyBypOverPt(false);
    analysiscuts->SetUsePreSelect(1);

    analysiscuts->SetUsePID(true);
    AliAODPidHF* PidHF = NULL;
    if(fUseStrongPID) {
        analysiscuts->SetUseStrongPid(3);
        analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
        analysiscuts->SetMaxPStrongPidK(1);
        analysiscuts->SetMaxPStrongPidpi(1);
    }

    analysiscuts->SetUseImpParProdCorrCut(false);
    analysiscuts->SetOptPileup(false);
    if(!fIsMC) {
        analysiscuts->SetMinCentrality(minc);
        analysiscuts->SetMaxCentrality(maxc);
    }

    analysiscuts->SetRemoveTrackletOutliers(true);//added on June 28
    analysiscuts->SetCutOnzVertexSPD(3);//needed for Pb-Pb 2015

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kCentral);
    if(fIsMC) analysiscuts->SetTriggerMask(AliVEvent::kMB);

    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb010);

    analysiscuts->SetMinPtCandidate(2.);
    analysiscuts->SetMaxPtCandidate(50.);

    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;

    TString triggername="kINT7_kCentral";
    if(fIsMC) triggername="kMB";
    TString PIDsuffix="";
    if(fUseStrongPID) PIDsuffix=Form("_strongPIDpt%0.f",maxPtstrongPID);

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename=Form("DplustoKpipiCuts_010_loose%s_Raa_%s.root",PIDsuffix.Data(),triggername.Data());
    TFile* fout=new TFile(filename.Data(),"RECREATE");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi* MakeFileForCutsDplus010_FiltTreeCreator2018(bool fUseStrongPID=false, double maxPtstrongPID=3.0, bool fIsMC=false){

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMaxChi2PerClusterTPC(2.5);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0060*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
    esdTrackCuts->SetPtRange(0.6,1.e10);

    float minc=0.,maxc=10.;
    const int nptbins=2;
    float* ptbins;
    ptbins=new float[nptbins+1];

    ptbins[0]=2.;
    ptbins[1]=5.;
    ptbins[2]=50.;

    const int nvars=14;
    float** anacutsval;
    anacutsval=new float*[nvars];
    for(int ic=0;ic<nvars;ic++){anacutsval[ic]=new float[nptbins];}

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

    for(int ipt=0;ipt<nptbins;ipt++){
        anacutsval[0][ipt]=0.2; //minv
        anacutsval[1][ipt]=0.6; //ptK
        anacutsval[2][ipt]=0.6; //ptPi
        anacutsval[3][ipt]=0.; //d0K
        anacutsval[4][ipt]=0.; //d0Pi
        anacutsval[5][ipt]=0.; //dist12
        anacutsval[8][ipt]=0.0; //pM
        anacutsval[10][ipt]=0.0; //sumd02
        anacutsval[11][ipt]=10000000000.; //dca
    }

    anacutsval[6][0]=0.04;
    anacutsval[7][0]=0.05;
    anacutsval[9][0]=0.99;
    anacutsval[12][0]=5.;
    anacutsval[13][0]=0.99;

    anacutsval[6][1]=0.06;
    anacutsval[7][1]=0.08;
    anacutsval[9][1]=0.97;
    anacutsval[12][1]=3.;
    anacutsval[13][1]=0.97;

    AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetMinNumTPCClsForPID(50);
    analysiscuts->SetScaleNormDLxyBypOverPt(false);
    analysiscuts->SetUsePreSelect(1);

    analysiscuts->SetUsePID(true);
    AliAODPidHF* PidHF = NULL;
    if(fUseStrongPID) {
        analysiscuts->SetUseStrongPid(3);
        analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
        analysiscuts->SetMaxPStrongPidK(1);
        analysiscuts->SetMaxPStrongPidpi(1);
    }

    analysiscuts->SetUseImpParProdCorrCut(false);
    analysiscuts->SetOptPileup(false);
    if(!fIsMC) {
        analysiscuts->SetMinCentrality(minc);
        analysiscuts->SetMaxCentrality(maxc);
    }

    analysiscuts->SetRemoveTrackletOutliers(true);//added on June 28
    analysiscuts->SetCutOnzVertexSPD(3);//needed for Pb-Pb 2015

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kCentral);
    if(fIsMC) analysiscuts->SetTriggerMask(AliVEvent::kMB);

    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb010);

    analysiscuts->SetMinPtCandidate(2.);
    analysiscuts->SetMaxPtCandidate(50.);

    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;

    TString triggername="kINT7_kCentral";
    if(fIsMC) triggername="kMB";
    TString PIDsuffix="";
    if(fUseStrongPID) PIDsuffix=Form("_strongPIDpt%0.f",maxPtstrongPID);

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename=Form("DplustoKpipiCuts_010_filttreecreator%s_Raa_%s.root",PIDsuffix.Data(),triggername.Data());
    TFile* fout=new TFile(filename.Data(),"RECREATE");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}
