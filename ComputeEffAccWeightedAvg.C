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
#include "TMath.h"
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


        for(int iPt=0; iPt<effC[0]->GetNbinsX(); iPt++) {
            double weightEffC = 0., weightUncEffC = 0., weightEffB = 0., weightUncEffB = 0., sumOfBR = 0.;
            for(int iCh=0; iCh<4; iCh++) {
                weightEffC += effC[iCh]->GetBinContent(iPt + 1) * BR[iCh];
                weightUncEffC += effC[iCh]->GetBinError(iPt + 1) * effC[iCh]->GetBinError(iPt + 1) * BR[iCh] * BR[iCh];
                weightEffB += effB[iCh]->GetBinContent(iPt + 1) * BR[iCh];
                weightUncEffB += effB[iCh]->GetBinError(iPt + 1) * effB[iCh]->GetBinError(iPt + 1) * BR[iCh] * BR[iCh];
                sumOfBR += BR[iCh];
            }
            weightEffC /= sumOfBR;
            weightUncEffC = TMath::Sqrt(weightUncEffC) / sumOfBR;
            weightEffB /= sumOfBR;
            weightUncEffB = TMath::Sqrt(weightUncEffB) / sumOfBR;
            effC[iCh]->SetBinContent(iPt + 1, weightEffC);
            effC[iCh]->SetBinError(iPt + 1, weightUncEffC);
            effB[iCh]->SetBinContent(iPt + 1, weightEffB);
            effB[iCh]->SetBinError(iPt + 1, weightUncEffB);
        }

        TFile outFile(Form("%s/%s", effdir.Data(),outFileName.Data()),"recreate");
        effCw->Write("hAccEffPrompt");
        effBw->Write("hAccEffFD");
    }

    return 0;

}