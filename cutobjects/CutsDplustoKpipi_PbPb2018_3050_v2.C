#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>
enum etaregion{kEtaPos,kEtaNeg,kEtaFull};

//__________________________________________________________________________________________
void MakeFileForCuts_Central(){

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                         AliESDtrackCuts::kAny);
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.4,1.e10);

  TString cent="";
  //centrality selection (Pb-Pb)
  Float_t minc=30.,maxc=50;
  const Int_t nptbins=15;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];

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

  const Int_t nvars=14;
  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}


  Int_t ic=0;//minv
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.2;
  }

  ic=1;//ptK
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[1][ipt]=0.4;
  }

  ic=2;//ptPi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.4;
  }

  ic=3;//d0K
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=4;//d0Pi
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=5;//dist12
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }

  ic=6;//sigvert
  anacutsval[ic][0]=0.020;//2.0-3.0
  anacutsval[ic][1]=0.022;//3.0-4.0
  anacutsval[ic][2]=0.022;//4.0-5.0
  anacutsval[ic][3]=0.022;//5.0-6.0
  anacutsval[ic][4]=0.024;//6.0-7.0
  anacutsval[ic][5]=0.024;//7.0-8.0
  anacutsval[ic][6]=0.024;//8.0-9.0
  anacutsval[ic][7]=0.024;//9.0-10.0
  anacutsval[ic][8]=0.024;//10.0-11.0
  anacutsval[ic][9]=0.024;//11.0-12.0
  anacutsval[ic][10]=0.024;//12.0-14.0
  anacutsval[ic][11]=0.024;//14.0-16.0
  anacutsval[ic][12]=0.024;//16.0-24.0
  anacutsval[ic][13]=0.034;//24.0-36.0
  anacutsval[ic][14]=0.034;//36.0-50.0

  ic=7;//declen
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.12;
  }
  anacutsval[ic][0]=0.08;//2.0-3.0
  anacutsval[ic][12]=0.16;//16.0-24.0
  anacutsval[ic][13]=0.16;//24.0-36.0
  anacutsval[ic][14]=0.16;//36.0-50.0

  ic=8;//pM
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }

  //cosp
  ic=9;
  anacutsval[ic][0]=0.996;//2.0-3.0
  anacutsval[ic][1]=0.996;//3.0-4.0
  anacutsval[ic][2]=0.995;//4.0-5.0
  anacutsval[ic][3]=0.995;//5.0-6.0
  anacutsval[ic][4]=0.995;//6.0-7.0
  anacutsval[ic][5]=0.995;//7.0-8.0
  anacutsval[ic][6]=0.990;//8.0-9.0
  anacutsval[ic][7]=0.990;//9.0-10.0
  anacutsval[ic][8]=0.990;//10.0-11.0
  anacutsval[ic][9]=0.990;//11.0-12.0
  anacutsval[ic][10]=0.990;//12.0-14.0
  anacutsval[ic][11]=0.990;//14.0-16.0
  anacutsval[ic][12]=0.980;//16.0-24.0
  anacutsval[ic][13]=0.970;//24.0-36.0
  anacutsval[ic][14]=0.950;//36.0-50.0

  ic=10;//sumd02
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.0;
  }

  ic=11;//dca
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=10000000000.;
  }

  ic=12;//ndlXY
  anacutsval[ic][0]=12.;//2.0-3.0
  anacutsval[ic][1]=12.;//3.0-4.0
  anacutsval[ic][2]=12.;//4.0-5.0
  anacutsval[ic][3]=12.;//5.0-6.0
  anacutsval[ic][4]=10.;//6.0-7.0
  anacutsval[ic][5]=10.;//7.0-8.0
  anacutsval[ic][6]=10.;//8.0-9.0
  anacutsval[ic][7]=10.;//9.0-10.0
  anacutsval[ic][8]=9.;//10.0-11.0
  anacutsval[ic][9]=9.;//11.0-12.0
  anacutsval[ic][10]=9.;//12.0-14.0
  anacutsval[ic][11]=9.;//14.0-16.0
  anacutsval[ic][12]=8.;//16.0-24.0
  anacutsval[ic][13]=8.;//24.0-36.0
  anacutsval[ic][14]=6.;//36.0-50.0

  ic=13;//cospXY
  anacutsval[ic][0]=0.996;//2.0-3.0
  anacutsval[ic][1]=0.996;//3.0-4.0
  anacutsval[ic][2]=0.995;//4.0-5.0
  anacutsval[ic][3]=0.995;//5.0-6.0
  anacutsval[ic][4]=0.995;//6.0-7.0
  anacutsval[ic][5]=0.995;//7.0-8.0
  anacutsval[ic][6]=0.990;//8.0-9.0
  anacutsval[ic][7]=0.990;//9.0-10.0
  anacutsval[ic][8]=0.990;//10.0-11.0
  anacutsval[ic][9]=0.990;//11.0-12.0
  anacutsval[ic][10]=0.990;//12.0-14.0
  anacutsval[ic][11]=0.990;//14.0-16.0
  anacutsval[ic][12]=0.980;//16.0-24.0
  anacutsval[ic][13]=0.970;//24.0-36.0
  anacutsval[ic][14]=0.950;//36.0-50.0

  Float_t *d0cutsval=new Float_t[nptbins];
  for(Int_t ipt=0;ipt<nptbins;ipt++){ //d0
    d0cutsval[ipt]=60;
  }
  d0cutsval[0]=80;//2.0-3.0

  Float_t *d0d0expcutsval=new Float_t[nptbins];
  for(Int_t ipt=0;ipt<nptbins;ipt++){ //d0d0exp
    d0d0expcutsval[ipt]=2.5;
  }
  d0d0expcutsval[0]=1.5;//2.0-3.0
  d0d0expcutsval[1]=1.5;//3.0-4.0
  d0d0expcutsval[2]=2.0;//4.0-5.0
  d0d0expcutsval[3]=2.0;//5.0-6.0
  d0d0expcutsval[14]=3.0;//36.0-50.0

  AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
  analysiscuts->SetName("AnalysisCuts");
  analysiscuts->SetTitle("Cuts for Dplus Analysis and CF");
  analysiscuts->SetPtBins(nptbins+1,ptbins);
  analysiscuts->SetCuts(nvars,nptbins,anacutsval);
  analysiscuts->Setd0Cut(nptbins,d0cutsval);
  analysiscuts->Setd0MeasMinusExpCut(nptbins,d0d0expcutsval);
  analysiscuts->AddTrackCuts(esdTrackCuts);
  analysiscuts->SetScaleNormDLxyBypOverPt(kFALSE);
  analysiscuts->SetUsePreSelect(1);

  cout<<"************** checking old PID (it should be FALSE by default - July 10)--> "<<analysiscuts->GetPidHF()->GetOldPid()<<endl;

  // analysiscuts->SetUsePID(kFALSE);
  analysiscuts->SetUsePID(kTRUE);
  analysiscuts->SetUseImpParProdCorrCut(kFALSE);
  analysiscuts->SetOptPileup(kFALSE);
  analysiscuts->SetMinCentrality(minc);
  analysiscuts->SetMaxCentrality(maxc);

  analysiscuts->SetRemoveTrackletOutliers(kTRUE);//added on June 28
  analysiscuts->SetCutOnzVertexSPD(3);//needed for Pb-Pb 2015

  cent=Form("%.0f%.0f",minc,maxc);
  analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
  analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
  analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
  analysiscuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kSemiCentral);

  analysiscuts->SetMinPtCandidate(2.);
  analysiscuts->SetMaxPtCandidate(50.);

  cout<<"This is the object I'm going to save:"<<nptbins<<endl;

  analysiscuts->PrintAll();
  analysiscuts->PrintTrigger();
  TString filename = "DplustoKpipiCuts_2018_3050_central_v2_kINT7_kSemiCentral.root";
  TFile* fout=new TFile(filename.Data(),"RECREATE");
  fout->cd();
  analysiscuts->Write();
  fout->Close();
}