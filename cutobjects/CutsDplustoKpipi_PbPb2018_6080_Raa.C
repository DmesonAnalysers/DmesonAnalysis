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
// 1) MakeFileForCutsDplus6080_Central2018 --> central cuts of 2015 and 2018 analyses
// 2) MakeFileForCutsDplus6080_Loose2018 --> loose cuts of 2018 analysis (sparse)
// 3) MakeFileForCutsDplus6080_FiltTreeCreator2018 --> filtering cuts of 2018 analysis (tree creator)
//____________________________________________________________________________________________________//

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi* MakeFileForCutsDplus6080_Central2018(bool fUseStrongPID=false, double maxPtstrongPID=3.0, bool fIsMC=false, int addTrackCut = 0){

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.4,1.e10);

    TString cent="";
    float minc=60.,maxc=80.;
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
        anacutsval[1][ipt]=0.2; //ptK
        anacutsval[2][ipt]=0.2; //ptPi
        anacutsval[3][ipt]=0.2; //d0K
        anacutsval[4][ipt]=0.2; //d0Pi
        anacutsval[5][ipt]=0.2; //dist12
        anacutsval[8][ipt]=0.0; //pM
        anacutsval[10][ipt]=0.0; //sumd02
        anacutsval[11][ipt]=10000000000.; //dca
    }

    int ic=6;//sigvert
    anacutsval[ic][0]=0.028;//2.0-3.0
    anacutsval[ic][1]=0.028;//3.0-4.0
    anacutsval[ic][2]=0.028;//4.0-5.0
    anacutsval[ic][3]=0.028;//5.0-6.0
    anacutsval[ic][4]=0.028;//6.0-7.0
    anacutsval[ic][5]=0.028;//7.0-8.0
    anacutsval[ic][6]=0.030;//8.0-9.0
    anacutsval[ic][7]=0.030;//9.0-10.0
    anacutsval[ic][8]=0.030;//10.0-11.0
    anacutsval[ic][9]=0.030;//11.0-12.0
    anacutsval[ic][10]=0.030;//12.0-14.0
    anacutsval[ic][11]=0.030;//14.0-16.0
    anacutsval[ic][12]=0.034;//16.0-24.0
    anacutsval[ic][13]=0.050;//24.0-36.0
    anacutsval[ic][14]=0.050;//36.0-50.0

    ic=7;//declen
    anacutsval[ic][0]=0.06;//2.0-3.0
    anacutsval[ic][1]=0.06;//3.0-4.0
    anacutsval[ic][2]=0.06;//4.0-5.0
    anacutsval[ic][3]=0.08;//5.0-6.0
    anacutsval[ic][4]=0.10;//6.0-7.0
    anacutsval[ic][5]=0.10;//7.0-8.0
    anacutsval[ic][6]=0.10;//8.0-9.0
    anacutsval[ic][7]=0.10;//9.0-10.0
    anacutsval[ic][8]=0.12;//10.0-11.0
    anacutsval[ic][9]=0.12;//11.0-12.0
    anacutsval[ic][10]=0.12;//12.0-14.0
    anacutsval[ic][11]=0.12;//14.0-16.0
    anacutsval[ic][12]=0.12;//16.0-24.0
    anacutsval[ic][13]=0.12;//24.0-36.0
    anacutsval[ic][14]=0.12;//36.0-50.0

    //cosp
    ic=9;
    anacutsval[ic][0]=0.992;//2.0-3.0
    anacutsval[ic][1]=0.992;//3.0-4.0
    anacutsval[ic][2]=0.992;//4.0-5.0
    anacutsval[ic][3]=0.992;//5.0-6.0
    anacutsval[ic][4]=0.990;//6.0-7.0
    anacutsval[ic][5]=0.990;//7.0-8.0
    anacutsval[ic][6]=0.985;//8.0-9.0
    anacutsval[ic][7]=0.985;//9.0-10.0
    anacutsval[ic][8]=0.985;//10.0-11.0
    anacutsval[ic][9]=0.985;//11.0-12.0
    anacutsval[ic][10]=0.980;//12.0-14.0
    anacutsval[ic][11]=0.980;//14.0-16.0
    anacutsval[ic][12]=0.970;//16.0-24.0
    anacutsval[ic][13]=0.970;//24.0-36.0
    anacutsval[ic][14]=0.970;//36.0-50.0

    ic=12;//ndlXY
    anacutsval[ic][0]=12.;//2.0-3.0
    anacutsval[ic][1]=12.;//3.0-4.0
    anacutsval[ic][2]=11.;//4.0-5.0
    anacutsval[ic][3]=10.;//5.0-6.0
    anacutsval[ic][4]=10.;//6.0-7.0
    anacutsval[ic][5]=10.;//7.0-8.0
    anacutsval[ic][6]=7.;//8.0-9.0
    anacutsval[ic][7]=7.;//9.0-10.0
    anacutsval[ic][8]=7.;//10.0-11.0
    anacutsval[ic][9]=7.;//11.0-12.0
    anacutsval[ic][10]=6.;//12.0-14.0
    anacutsval[ic][11]=6.;//14.0-16.0
    anacutsval[ic][12]=4.;//16.0-24.0
    anacutsval[ic][13]=4.;//24.0-36.0
    anacutsval[ic][14]=4.;//36.0-50.0

    ic=13;//cospXY
    anacutsval[ic][0]=0.994;//2.0-3.0
    anacutsval[ic][1]=0.994;//3.0-4.0
    anacutsval[ic][2]=0.992;//4.0-5.0
    anacutsval[ic][3]=0.992;//5.0-6.0
    anacutsval[ic][4]=0.992;//6.0-7.0
    anacutsval[ic][5]=0.992;//7.0-8.0
    anacutsval[ic][6]=0.990;//8.0-9.0
    anacutsval[ic][7]=0.990;//9.0-10.0
    anacutsval[ic][8]=0.990;//10.0-11.0
    anacutsval[ic][9]=0.990;//11.0-12.0
    anacutsval[ic][10]=0.985;//12.0-14.0
    anacutsval[ic][11]=0.985;//14.0-16.0
    anacutsval[ic][12]=0.970;//16.0-24.0
    anacutsval[ic][13]=0.970;//24.0-36.0
    anacutsval[ic][14]=0.970;//36.0-50.0

    Float_t *d0cutsval=new Float_t[nptbins];
    for(Int_t ipt=0;ipt<nptbins;ipt++){ //d0
      d0cutsval[ipt]=80;
    }
    d0cutsval[0]=100;//2.0-3.0
    d0cutsval[1]=100;//3.0-4.0
    d0cutsval[2]=100;//4.0-5.0
    d0cutsval[3]=100;//5.0-6.0
    d0cutsval[4]=90;//6.0-7.0

    Float_t *d0d0expcutsval=new Float_t[nptbins];
    for(Int_t ipt=0;ipt<nptbins;ipt++){ //d0d0exp
      d0d0expcutsval[ipt]=2.5;
    }
    d0d0expcutsval[0]=2.0;//2.0-3.0
    d0d0expcutsval[1]=2.0;//3.0-4.0
    d0d0expcutsval[10]=3.0;//12.0-14.0
    d0d0expcutsval[11]=3.0;//14.0-16.0
    d0d0expcutsval[12]=3.0;//16.0-24.0
    d0d0expcutsval[13]=3.0;//24.0-36.0
    d0d0expcutsval[14]=3.0;//36.0-50.0

    AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->Setd0Cut(nptbins,d0cutsval);
    analysiscuts->Setd0MeasMinusExpCut(nptbins,d0d0expcutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetScaleNormDLxyBypOverPt(false);

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

    cent=Form("%.0f%.0f",minc,maxc);
    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    if(fIsMC) analysiscuts->SetTriggerMask(AliVEvent::kMB);
    
    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb6080);

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

    analysiscuts->SetMinPtCandidate(2.);
    analysiscuts->SetMaxPtCandidate(50.);

    cout<<"This is the object I'm going to save:"<<nptbins<<endl;

    TString triggername="kINT7";
    if(fIsMC) triggername="kMB";
    TString PIDsuffix="";
    if(fUseStrongPID) PIDsuffix=Form("_strongPIDpt%0.f",maxPtstrongPID);

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename=Form("DplustoKpipiCuts_6080_central%s_Raa_%s%s.root",PIDsuffix.Data(),triggername.Data(),trackCutName.Data());
    TFile* fout=new TFile(filename.Data(),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi* MakeFileForCutsDplus6080_Loose2018(bool fUseStrongPID=false, double maxPtstrongPID=3.0, bool fIsMC=false){

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.4,1.e10);

    TString cent="";
    float minc=60.,maxc=80.;
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
        anacutsval[1][ipt]=0.2; //ptK
        anacutsval[2][ipt]=0.2; //ptPi
        anacutsval[3][ipt]=0.2; //d0K
        anacutsval[4][ipt]=0.2; //d0Pi
        anacutsval[5][ipt]=0.2; //dist12
        anacutsval[8][ipt]=0.0; //pM
        anacutsval[10][ipt]=0.0; //sumd02
        anacutsval[11][ipt]=10000000000.; //dca
    }

    int ic=6;//sigvert
    anacutsval[ic][0]=0.040;//2.0-3.0
    anacutsval[ic][1]=0.040;//3.0-4.0
    anacutsval[ic][2]=0.040;//4.0-5.0
    anacutsval[ic][3]=0.040;//5.0-6.0
    anacutsval[ic][4]=0.040;//6.0-7.0
    anacutsval[ic][5]=0.040;//7.0-8.0
    anacutsval[ic][6]=0.040;//8.0-9.0
    anacutsval[ic][7]=0.040;//9.0-10.0
    anacutsval[ic][8]=0.040;//10.0-11.0
    anacutsval[ic][9]=0.040;//11.0-12.0
    anacutsval[ic][10]=0.050;//12.0-14.0
    anacutsval[ic][11]=0.050;//14.0-16.0
    anacutsval[ic][12]=0.050;//16.0-24.0
    anacutsval[ic][13]=0.060;//24.0-36.0
    anacutsval[ic][14]=0.060;//36.0-50.0

    ic=7;//declen
    for(Int_t ipt=0;ipt<nptbins;ipt++){
        anacutsval[ic][ipt]=0.04;
    }

    //cosp
    ic=9;
    anacutsval[ic][0]=0.980;//2.0-3.0
    anacutsval[ic][1]=0.980;//3.0-4.0
    anacutsval[ic][2]=0.980;//4.0-5.0
    anacutsval[ic][3]=0.980;//5.0-6.0
    anacutsval[ic][4]=0.980;//6.0-7.0
    anacutsval[ic][5]=0.980;//7.0-8.0
    anacutsval[ic][6]=0.970;//8.0-9.0
    anacutsval[ic][7]=0.970;//9.0-10.0
    anacutsval[ic][8]=0.970;//10.0-11.0
    anacutsval[ic][9]=0.970;//11.0-12.0
    anacutsval[ic][10]=0.950;//12.0-14.0
    anacutsval[ic][11]=0.950;//14.0-16.0
    anacutsval[ic][12]=0.950;//16.0-24.0
    anacutsval[ic][13]=0.950;//24.0-36.0
    anacutsval[ic][14]=0.950;//36.0-50.0

    ic=12;//ndlXY
    anacutsval[ic][0]=5.;//2.0-3.0
    anacutsval[ic][1]=5.;//3.0-4.0
    anacutsval[ic][2]=5.;//4.0-5.0
    anacutsval[ic][3]=5.;//5.0-6.0
    anacutsval[ic][4]=5.;//6.0-7.0
    anacutsval[ic][5]=5.;//7.0-8.0
    anacutsval[ic][6]=5.;//8.0-9.0
    anacutsval[ic][7]=4.;//9.0-10.0
    anacutsval[ic][8]=4.;//10.0-11.0
    anacutsval[ic][9]=4.;//11.0-12.0
    anacutsval[ic][10]=4.;//12.0-14.0
    anacutsval[ic][11]=4.;//14.0-16.0
    anacutsval[ic][12]=4.;//16.0-24.0
    anacutsval[ic][13]=4.;//24.0-36.0
    anacutsval[ic][14]=4.;//36.0-50.0

    ic=13;//cospXY
    anacutsval[ic][0]=0.980;//2.0-3.0
    anacutsval[ic][1]=0.980;//3.0-4.0
    anacutsval[ic][2]=0.980;//4.0-5.0
    anacutsval[ic][3]=0.980;//5.0-6.0
    anacutsval[ic][4]=0.980;//6.0-7.0
    anacutsval[ic][5]=0.980;//7.0-8.0
    anacutsval[ic][6]=0.970;//8.0-9.0
    anacutsval[ic][7]=0.970;//9.0-10.0
    anacutsval[ic][8]=0.970;//10.0-11.0
    anacutsval[ic][9]=0.970;//11.0-12.0
    anacutsval[ic][10]=0.950;//12.0-14.0
    anacutsval[ic][11]=0.950;//14.0-16.0
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

    cent=Form("%.0f%.0f",minc,maxc);
    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    if(fIsMC) analysiscuts->SetTriggerMask(AliVEvent::kMB);
    
    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb6080);

    analysiscuts->SetMinPtCandidate(2.);
    analysiscuts->SetMaxPtCandidate(50.);

    cout<<"This is the object I'm going to save:"<<nptbins<<endl;

    TString triggername="kINT7";
    if(fIsMC) triggername="kMB";
    TString PIDsuffix="";
    if(fUseStrongPID) PIDsuffix=Form("_strongPIDpt%0.f",maxPtstrongPID);

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename=Form("DplustoKpipiCuts_6080_loose%s_Raa_%s.root",PIDsuffix.Data(),triggername.Data());
    TFile* fout=new TFile(filename.Data(),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}

//__________________________________________________________________________________________
AliRDHFCutsDplustoKpipi* MakeFileForCutsDplus6080_FiltTreeCreator2018(bool fUseStrongPID=false, double maxPtstrongPID=3.0, bool fIsMC=false){

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false);
    esdTrackCuts->SetRequireTPCRefit(true);
    esdTrackCuts->SetRequireITSRefit(true);
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.4,1.e10);

    TString cent="";
    float minc=60.,maxc=80.;
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
        anacutsval[1][ipt]=0.2; //ptK
        anacutsval[2][ipt]=0.2; //ptPi
        anacutsval[3][ipt]=0.2; //d0K
        anacutsval[4][ipt]=0.2; //d0Pi
        anacutsval[5][ipt]=0.2; //dist12
        anacutsval[8][ipt]=0.0; //pM
        anacutsval[10][ipt]=0.0; //sumd02
        anacutsval[11][ipt]=10000000000.; //dca
    }

    anacutsval[6][0]=0.03;
    anacutsval[7][0]=0.03;
    anacutsval[9][0]=0.98;
    anacutsval[12][0]=3.;
    anacutsval[13][0]=0.98;

    anacutsval[6][0]=0.03;
    anacutsval[7][0]=0.08;
    anacutsval[9][0]=0.95;
    anacutsval[12][0]=2.;
    anacutsval[13][0]=0.95;

    AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
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

    cent=Form("%.0f%.0f",minc,maxc);
    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    if(fIsMC) analysiscuts->SetTriggerMask(AliVEvent::kMB);
    
    if(!fIsMC)
        analysiscuts->SetUseTimeRangeCutForPbPb2018(true);
    if(!fIsMC)
      analysiscuts->EnableNsigmaDataDrivenCorrection(true,AliAODPidHF::kPbPb6080);

    analysiscuts->SetMinPtCandidate(2.);
    analysiscuts->SetMaxPtCandidate(50.);

    cout<<"This is the object I'm going to save:"<<nptbins<<endl;

    TString triggername="kINT7";
    if(fIsMC) triggername="kMB";
    TString PIDsuffix="";
    if(fUseStrongPID) PIDsuffix=Form("_strongPIDpt%0.f",maxPtstrongPID);

    analysiscuts->PrintAll();
    analysiscuts->PrintTrigger();
    TString filename=Form("DplustoKpipiCuts_6080_filttreecreator%s_Raa_%s.root",PIDsuffix.Data(),triggername.Data());
    TFile* fout=new TFile(filename.Data(),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    return analysiscuts;
}
