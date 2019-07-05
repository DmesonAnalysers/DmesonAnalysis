#include <iostream>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>
#include <AliRDHFCutsDstoKKpi.h>

void MakeFileForCuts_Loose2018(Bool_t fUseStrongPID = kTRUE, Double_t maxPtstrongPID = 8.0, Bool_t fIsMC=kFALSE) {

    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    //default
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    //esdTrackCuts->SetMinNClustersITS(4); // default is 5
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny);
    // default is kBoth, otherwise kAny
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3,1.e10);
    
    Float_t mincen=60.;
    Float_t maxcen=80.;

    if(fIsMC){
      mincen=0.;
      maxcen=100.;
    }
    
    const Int_t nptbins=9;
    Float_t* ptbins;
    ptbins=new Float_t[nptbins+1];
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

    const Int_t nvars=20;
    
    Float_t** anacutsval;
    anacutsval=new Float_t*[nvars];
    
    for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
    for(Int_t ipt=0;ipt<nptbins;ipt++){
        
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
    anacutsval[7][1]=0.03;   //decay length
    anacutsval[9][1]=0.96;   //cosP
    anacutsval[12][1]=0.010; //Mass Phi
    anacutsval[14][1]=0.0;  //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][1]=1.0;  //CosP labFrame
    anacutsval[16][1]=0.;  //decayXY
    anacutsval[17][1]=0.;    //normdecay
    anacutsval[18][1]=3.0;   //normdecayXY
    anacutsval[19][1]=0.96;  //CosPXY
    
    anacutsval[6][2]=0.045;   //sigmavert
    anacutsval[7][2]=0.03;   //decay length
    anacutsval[9][2]=0.96;   //cosP
    anacutsval[12][2]=0.008; //Mass Phi
    anacutsval[14][2]=0.0;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][2]=1.0;   //CosP labFrame
    anacutsval[16][2]=0.;  //decayXY
    anacutsval[17][2]=0.;    //normdecay
    anacutsval[18][2]=3.0;   //normdecayXY
    anacutsval[19][2]=0.96;  //CosPXY
    
    anacutsval[6][3]=0.045;   //sigmavert
    anacutsval[7][3]=0.03;   //decay length
    anacutsval[9][3]=0.96;   //cosP
    anacutsval[12][3]=0.010; //Mass Phi
    anacutsval[14][3]=0.0;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][3]=1.0;  //CosP labFrame
    anacutsval[16][3]=0.;   //decayXY
    anacutsval[17][3]=0.;    //normdecay
    anacutsval[18][3]=3.0;   //normdecayXY
    anacutsval[19][3]=0.96;  //CosPXY

    anacutsval[6][4]=0.045;   //sigmavert
    anacutsval[7][4]=0.03;   //decay length
    anacutsval[9][4]=0.95;   //cosP
    anacutsval[12][4]=0.010; //Mass Phi
    anacutsval[14][4]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][4]=1.0;  //CosP labFrame
    anacutsval[16][4]=0.;   //decayXY
    anacutsval[17][4]=0.;    //normdecay
    anacutsval[18][4]=3.0;   //normdecayXY
    anacutsval[19][4]=0.95;  //CosPXY

    anacutsval[6][5]=0.045;   //sigmavert
    anacutsval[7][5]=0.04;   //decay length
    anacutsval[9][5]=0.95;   //cosP
    anacutsval[12][5]=0.010; //Mass Phi
    anacutsval[14][5]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][5]=1.0;  //CosP labFrame
    anacutsval[16][5]=0.;   //decayXY
    anacutsval[17][5]=0.;    //normdecay
    anacutsval[18][5]=3.0;   //normdecayXY
    anacutsval[19][5]=0.95;  //CosPXY

    anacutsval[6][6]=0.045;   //sigmavert
    anacutsval[7][6]=0.04;   //decay length
    anacutsval[9][6]=0.95;   //cosP
    anacutsval[12][6]=0.010; //Mass Phi
    anacutsval[14][6]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][6]=1.0;  //CosP labFrame
    anacutsval[16][6]=0.;   //decayXY
    anacutsval[17][6]=0.;    //normdecay
    anacutsval[18][6]=3.0;   //normdecayXY
    anacutsval[19][6]=0.95;  //CosPXY

    anacutsval[6][7]=0.055;   //sigmavert
    anacutsval[7][7]=0.03;   //decay length
    anacutsval[9][7]=0.95;   //cosP
    anacutsval[12][7]=0.010; //Mass Phi
    anacutsval[14][7]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][7]=1.0;  //CosP labFrame
    anacutsval[16][7]=0.;   //decayXY
    anacutsval[17][7]=0.;    //normdecay
    anacutsval[18][7]=3.0;   //normdecayXY
    anacutsval[19][7]=0.95;  //CosPXY

    anacutsval[6][8]=0.055;   //sigmavert
    anacutsval[7][8]=0.03;   //decay length
    anacutsval[9][8]=0.95;   //cosP
    anacutsval[12][8]=0.010; //Mass Phi
    anacutsval[14][8]=0.;   //Abs(CosineKpiPhiRFrame)^3
    anacutsval[15][8]=1.0;  //CosP labFrame
    anacutsval[16][8]=0.;   //decayXY
    anacutsval[17][8]=0.;    //normdecay
    anacutsval[18][8]=3.0;   //normdecayXY
    anacutsval[19][8]=0.95;  //CosPXY

    Float_t topomCuts[nptbins] = {4.,4.,4.,5.,5.,5.,6.,6.,6.}; //topomatic
    Float_t d0Cuts[nptbins] = {500.,500.,500.,500.,500.,500.,500.,500.,500.}; //impparXY
    
    AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Ds Analysis and CF");
    analysiscuts->SetUsePreSelect(1);
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->Setd0MeasMinusExpCut(nptbins,topomCuts);
    analysiscuts->Setd0Cut(nptbins,d0Cuts);
    
    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    if(!fIsMC) 
      analysiscuts->SetTriggerMask(AliVEvent::kINT7);
    else 
      analysiscuts->SetTriggerMask(AliVEvent::kMB);
    
    analysiscuts->SetUsePID(kTRUE);
    if(fUseStrongPID) {
      analysiscuts->SetPidOption(1); //0=kConservative,1=kStrong
      analysiscuts->SetMaxPtStrongPid(maxPtstrongPID);
    }
    else analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong
    // if(!fIsMC)
    //   analysiscuts->EnableNsigmaDataDrivenCorrection(kTRUE,AliAODPidHF::kPbPb3050);

    analysiscuts->SetOptPileup(kFALSE);
    analysiscuts->SetMinCentrality(mincen);
    analysiscuts->SetMaxCentrality(maxcen);
    
    analysiscuts->SetMinPtCandidate(2.);
    analysiscuts->SetMaxPtCandidate(36.);
    analysiscuts->SetRemoveDaughtersFromPrim(kFALSE);
  
    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;
    
    analysiscuts->PrintAll();
    TString triggername = "kINT7";
    if(fIsMC)
     triggername = "kMB";
    TFile* fout=new TFile(Form("DstoKKpiCuts_6080_loose_strongPIDpt%0.f_Raa_%s.root",maxPtstrongPID, triggername.Data()),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();
}