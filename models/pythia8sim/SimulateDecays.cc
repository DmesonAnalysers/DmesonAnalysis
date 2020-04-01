//_________________________________________________________//
// Macro for the simulation of B hadron decays             //
// Main Function: SimulateDecays                           //
//_________________________________________________________//

#if !defined (__CINT__) || defined (__CLING__)

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>

#include "yaml-cpp/yaml.h"

#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TTree.h>
#include <TH1F.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TSystem.h>
#include <TStyle.h>
#include <AliDecayerPythia8.h>
#include <TPaveStats.h>

#endif

using namespace std;

enum predOpt{kCentral, kMin, kMax};

//____________________________________________________________________________________
//function prototypes
void SimulateDecays(TString cfgFileName = "config_Bhadron_decay.yml");
vector<TParticle*> SearchForDaughters(TClonesArray* array, vector<int> pdgCodeD);
TH1F* ReadFONLLPtShapeFile(string txtFileName = "", int whichPred=kCentral, double ptMin=0.025, double binWidth=0.05, int nBins=1001);

//____________________________________________________________________________________
//function implementations
void SimulateDecays(TString cfgFileName)
{
    //upload config
    YAML::Node config = YAML::LoadFile(cfgFileName.Data());
    vector<int> pdgCodeB = config["pdgCodeB"].as<vector<int> >();
    vector<int> pdgCodeD = config["pdgCodeD"].as<vector<int> >();
    int nDecaysPerSpecies = config["nDecaysPerSpecies"].as<int>();
    string outFileName = config["outFileName"].as<string>();
    string ptShapeBmesonFile = "";
    double ptMin = 0., ptBinWidth = -1.;
    int nPtBins = -1;
    int predOpt = -1;
    TH1F* hPtBShape = nullptr;

    cout << static_cast<bool>(config["BhadronPtShape"]["activate"].as<int>()) << endl;
    if(static_cast<bool>(config["BhadronPtShape"]["activate"].as<int>()))
    {
        ptShapeBmesonFile = config["BhadronPtShape"]["inFile"].as<string>();
        ptMin = config["BhadronPtShape"]["ptMin"].as<double>();
        ptBinWidth = config["BhadronPtShape"]["ptBinWidth"].as<double>();
        nPtBins = config["BhadronPtShape"]["nPtBins"].as<int>();
        string sPredOpt = config["BhadronPtShape"]["predOpt"].as<string>();
        if(sPredOpt == "kCentral")
            predOpt = kCentral;
        else if(sPredOpt == "kMin")
            predOpt = kMin;
        else if(sPredOpt == "kMax")
            predOpt = kMax;

        hPtBShape = ReadFONLLPtShapeFile(ptShapeBmesonFile, predOpt, ptMin, ptBinWidth, nPtBins);
    }

    //define decayer
    AliDecayerPythia8 *pdec = new AliDecayerPythia8();
    pdec->Init();

    //define tree
    TTree* fTreeDecays = new TTree("fTreeDecays", "fTreeDecays");

    int pdgB = -9999;
    float ptB = -1.;
    float pB = -1.;
    float yB = -1.;
    vector<float> ptD;
    vector<float> pD;
    vector<float> yD;
    vector<int> pdgD;
    float norm = 1.;
    if(static_cast<bool>(config["BhadronPtShape"]["activate"].as<int>()))
        norm = hPtBShape->Integral();

    fTreeDecays->Branch("pdgB", &pdgB);
    fTreeDecays->Branch("ptB", &ptB);
    fTreeDecays->Branch("pB", &pB);
    fTreeDecays->Branch("yB", &yB);
    fTreeDecays->Branch("ptD", &ptD);
    fTreeDecays->Branch("pD", &pD);
    fTreeDecays->Branch("yD", &yD);
    fTreeDecays->Branch("pdgD", &pdgD);
    fTreeDecays->Branch("norm", &norm);

    //decay simulation
    TRandom3 *gener = new TRandom3(0);
    TLorentzVector *vec = new TLorentzVector();
    TClonesArray *array = new TClonesArray("TParticle", 100);
    TDatabasePDG* db = TDatabasePDG::Instance();

    for (auto pdgCode: pdgCodeB)
    {
        int iDecay = 0;
        while (iDecay < nDecaysPerSpecies)
        {
            float massB = db->GetParticle(pdgCode)->Mass();
            if(!static_cast<bool>(config["BhadronPtShape"]["activate"].as<int>())) //decay at rest
                vec->SetPxPyPzE(0, 0, 0, massB); 
            else
            {
                ptB = hPtBShape->GetRandom();
                double phiB = gRandom->Rndm()*2*TMath::Pi();
                double yB = gRandom->Rndm()*2.-1.; // flat in -1<y<1
                double pxB = ptB * TMath::Cos(phiB);
                double pyB = ptB * TMath::Sin(phiB);
                double mtB = TMath::Sqrt(massB * massB + ptB * ptB);
                double pzB = mtB * TMath::SinH(yB);
                double EB = TMath::Sqrt(massB * massB + pxB * pxB + pyB * pyB + pzB * pzB);
                vec->SetPxPyPzE(pxB, pyB, pzB, EB); 
            }
            pdec->Decay(pdgCode, vec);
            array->Clear();

            int nentries = pdec->ImportParticles(array);
            TParticle *Bmes = dynamic_cast<TParticle *>(array->At(0));
            pdgB = Bmes->GetPdgCode();
            ptB = Bmes->Pt();
            pB = Bmes->P();
            yB = Bmes->Y();
            vector<TParticle*> dautokeep = SearchForDaughters(array, pdgCodeD);
            for(unsigned int iDau=0; iDau<dautokeep.size(); iDau++)
            {
                int pdgCodeDau = dautokeep[iDau]->GetPdgCode();
                ptD.push_back(dautokeep[iDau]->Pt());
                pD.push_back(dautokeep[iDau]->P());
                yD.push_back(dautokeep[iDau]->Y());
                pdgD.push_back(pdgCodeDau);
            }
            if(ptD.size() == 0)
            {
                ptD.push_back(-1);
                pD.push_back(-1);
                yD.push_back(-1);
                pdgD.push_back(-1);
            }
            else
            {
                if (iDecay % 100000 == 0)
                    cout << "Desired generation number " << iDecay << " for hadron " << pdgB << endl;
                iDecay++;
            }

            fTreeDecays->Fill();
            pdgB = -1;
            ptD.clear();
            pD.clear();
            yD.clear();
            pdgD.clear();
        }
    }

    TFile outfile(outFileName.data(), "recreate");
    fTreeDecays->Write();
    outfile.Close();
}

//________________________________________________________________________
vector<TParticle*> SearchForDaughters(TClonesArray* array, vector<int> pdgCodeD)
{
    vector<TParticle*> parttokeep;
    for (int iDau = 1; iDau < array->GetEntriesFast(); iDau++) //the first one is the B
    {
        TParticle* dau = dynamic_cast<TParticle *>(array->At(iDau));
        if(dau == nullptr)
            continue;
        int pdgDau = TMath::Abs(dau->GetPdgCode());
        if(find(pdgCodeD.begin(), pdgCodeD.end(), pdgDau) != pdgCodeD.end())
            parttokeep.push_back(dau);
    }

    return parttokeep;
}

//________________________________________________________________________
TH1F* ReadFONLLPtShapeFile(string txtFileName, int whichPred, double ptMin, double binWidth, int nBins)
{
    if(txtFileName.find("txt") == string::npos && txtFileName.find("dat") == string::npos && txtFileName.find("csv") == string::npos) {
        cerr << "ERROR: Wrong file format! Exit." << endl;
        return nullptr;
    }

    ifstream inSet(txtFileName.data());
    if(!inSet) {
        cerr << "ERROR: Please check if "<< txtFileName.data() << " is the right path. Exit." << endl;
        return nullptr;
    }

    vector<string> values;
    
    TH1F* hPtBShape = new TH1F("hPtBShape", "", nBins, ptMin, ptMin + nBins*binWidth);
    double centY = -1., minY = -1., maxY = -1.;
    int iPt = 1;
    while(!inSet.eof())
    {    
        vector<string> values;
        string line;
        getline(inSet, line);
        if(line.find("#") != string::npos || line.find("pt") != string::npos)
            continue;

        size_t pos = 0;
        while((pos = line.find(' ')) != string::npos)
        {
            values.push_back(line.substr(0, pos));
            if(values.size()==2) {
                stringstream convert(values[1]);
                if( !(convert >> centY) )
                    centY = -1;
            }
            else if(values.size()==3) {
                stringstream convert(values[2]);
                if( !(convert >> minY) )
                    minY = -1;
            }
            else if(values.size()==4) {
                stringstream convert(values[3]);
                if( !(convert >> maxY) )
                    maxY = -1;
            }
            line = line.substr(pos + 1);
        }
        if(whichPred == kCentral)
            hPtBShape->SetBinContent(iPt, centY);
        else if(whichPred == kMin)
            hPtBShape->SetBinContent(iPt, minY);
        else if(whichPred == kMax)
            hPtBShape->SetBinContent(iPt, maxY);
        iPt++;
    }

    inSet.close();

    return hPtBShape;
}