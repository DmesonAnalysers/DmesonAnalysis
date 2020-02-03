#if !defined (__CINT__) || defined (__CLING__)

#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>

#include "AliESDtrackCuts.h"
#include "AliRDHFCutsDstoKKpi.h"

#endif

//____________________________________________________________________________________________________//
// Methods:
// 1) MakeFileForCutsDs_FiltTree --> filtering cuts for tree used for ITS3 studies with ML
//____________________________________________________________________________________________________//


AliRDHFCutsDstoKKpi* MakeFileForCutsDs_FiltTree(double ptmin=1., double ptmax=50.) {


    AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(false); //cut track if sigma from track-to-vertex could not be calculated
    esdTrackCuts->SetRequireTPCRefit(true); // require TPC refit...what does refit means? 
    esdTrackCuts->SetRequireITSRefit(true); // require ITS refit...what does refit means?
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); // ??
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny); // requirement over the size and the form of the clusters for all the layers in ITS
    esdTrackCuts->SetMinDCAToVertexXY(0.);  //DCA should be the Detector Control Agent and should be able to control the status of the detector
    esdTrackCuts->SetPtRange(0.4,1.e10);    


    const int nptbins=3; 
    float ptbins[nptbins+1]={0.,5.,10.,50.};

    const int nvars = 20;

    float** anacutsval=new float*[nvars];
    for(int ic=0;ic<nvars;ic++){anacutsval[ic]=new float[nptbins];}

    /*
        Cut list
        0       "inv. mass [GeV]",
        1	    "pTK [GeV/c]",
        2	    "pTPi [GeV/c]",
        3	    "d0K [cm]",
        4	    "d0Pi [cm]",
        5		"dist12 [cm]",
        6		"sigmavert [cm]",
        7		"decLen [cm]",
        8		"ptMax [GeV/c]",
        9		"cosThetaPoint",
        10		"Sum d0^2 (cm^2)",
        11		"dca [cm]",
        12		"inv. mass (Mphi-MKK) [GeV]",
        13		"inv. mass (MKo*-MKpi) [GeV]",
        14   	"Abs(CosineKpiPhiRFrame)^3",
        15 		"CosPiDsLabFrame",
        16 		"DecLengthXY",
        17 		"NormDecayLength",
        18 		"NormDecayLengthXY",
        19 		"cosThetaPointXY"
     */
    
    anacutsval[7][0]=0.02;
    anacutsval[9][0]=0.97;
    anacutsval[12][0]=0.015;
    anacutsval[18][0]=4.;
    
    anacutsval[7][1]=0.02;
    anacutsval[9][1]=0.95;
    anacutsval[12][1]=0.015;
    anacutsval[18][1]=3.;
    
    anacutsval[12][2]=0.015;

    AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Ds Analysis and CF");
    analysiscuts->SetUsePreSelect(1);
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);

    analysiscuts->SetUseCentrality(AliRDHFCuts::kCentOff); //
    analysiscuts->SetTriggerClass("");//dont use for ppMB/ppMB_MC
    analysiscuts->ResetMaskAndEnableMBTrigger();//dont use for ppMB/ppMB_MC
    analysiscuts->SetTriggerMask(AliVEvent::kAny);

    analysiscuts->SetUsePID(true);
    analysiscuts->SetPidOption(0); //0=kConservative,1=kStrong

    analysiscuts->SetOptPileup(false);
    analysiscuts->SetMinPtCandidate(ptmin);
    analysiscuts->SetMaxPtCandidate(ptmax);
    analysiscuts->SetRemoveDaughtersFromPrim(false);

    std::cout<<"This is the object I'm going to save:"<<nptbins<<std::endl;

    analysiscuts->PrintAll();

    TFile* fout=new TFile(Form("DstoKKpiCutsITS3_filttree_kAny_consPID_pt%0.f_%0.f.root",ptmin,ptmax),"recreate");
    fout->cd();
    analysiscuts->Write();
    fout->Close();

    std::cout<<"\nCpu Time "<<stopwatch.CpuTime()<<"s\n"<<std::endl;

    return analysiscuts;
}
