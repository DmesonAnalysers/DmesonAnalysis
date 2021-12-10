//___________________________________________________________________________________//
// Macro for computing the weighted average of Lc resonant channels                  //
// Main Function: ComputeEffAccWeightedAvg                                           //
//___________________________________________________________________________________//

#if !defined (__CINT__) || defined (__CLING__)

#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include "TROOT.h"
#include "Riostream.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TDatabasePDG.h"

#include "AliHFInvMassFitter.h"
#include "AliVertexingHFUtils.h"

#endif

using namespace std;

int ComputeEffAccWeightedAvg(TString effdir, TString particle, TString cutset, TString outFileName){

    if(particle=="LctopKpi"){
        TString f[5] = {
            Form("%s/Eff_times_Acc_%s_NonRes%s.root", effdir.Data(), particle.Data(), cutset.Data()), //NonResonant
            Form("%s/Eff_times_Acc_%s_KStar%s.root", effdir.Data(), particle.Data(), cutset.Data()), //KStar
            Form("%s/Eff_times_Acc_%s_Delta%s.root", effdir.Data(), particle.Data(), cutset.Data()), //Delta
            Form("%s/Eff_times_Acc_%s_Lambda1520%s.root", effdir.Data(), particle.Data(), cutset.Data()), //Lambda1520
            Form("%s/Eff_times_Acc_%s%s.root", effdir.Data(), particle.Data(), cutset.Data()) //All
        };
        
        TFile *file[5];
        for (int i = 0; i < 5; i++){
            file[i] = TFile::Open(f[i].Data());
        }
        
        Double_t BR[4] = {
            //6.28 * 1e-02,        //all
            3.5 * 1e-02,         //non res
            1.96 * 0.667 * 1e-02, //Kstar
            1.08 * 1e-02,        //Delta
            2.2 * 0.225 * 1e-02 //L1520
        };

        TH1D *effC[5];
        TH1D *effB[5];
        for(int i = 0; i < 5; i++){
            effC[i] = (TH1D *)file[i]->Get("hAccEffPrompt");
            effB[i] = (TH1D *)file[i]->Get("hAccEffFD");
            cout << " n bins " << effC[i]->GetNbinsX() << endl;
        }
        TH1D *effCw = (TH1D *)effC[0]->Clone("hAccEffPrompt");
        TH1D *effBw = (TH1D *)effB[0]->Clone("hAccEffFD");
        effCw->Reset("icse"); 
        effBw->Reset("icse");


        for(int k=0; k<effC[0]->GetNbinsX(); k++){
            effCw->SetBinContent(k + 1, ((effC[0]->GetBinContent(k + 1) * BR[0]) + (effC[1]->GetBinContent(k + 1) * BR[1]) + (effC[2]->GetBinContent(k + 1) * BR[2]) + (effC[3]->GetBinContent(k + 1) * BR[3])) / (BR[0] + BR[1] + BR[2] + BR[3]));
            effCw->SetBinError(k + 1, effC[4]->GetBinError(k + 1));

            effBw->SetBinContent(k + 1, ((effB[0]->GetBinContent(k + 1) * BR[0]) + (effB[1]->GetBinContent(k + 1) * BR[1]) + (effB[2]->GetBinContent(k + 1) * BR[2]) + (effB[3]->GetBinContent(k + 1) * BR[3])) / (BR[0] + BR[1] + BR[2] + BR[3]));
            effBw->SetBinError(k + 1, effB[4]->GetBinError(k + 1));
        }

        TFile outFile(Form("%s/%s", effdir.Data(),outFileName.Data()),"recreate");
        effCw->Write("hAccEffPrompt");
        effBw->Write("hAccEffFD");
    }

    return 0;

}