#if !defined (__CINT__) || defined (__CLING__)

#include <vector>
#include <string>
#include <iostream>
#include "yaml-cpp/yaml.h"

#include <TFile.h>
#include <TTree.h>
#include "ROOT/RDataFrame.hxx"

#endif

using std::vector;
using std::string;
using std::cout;
using std::endl;

//______________________________________________________________________________________________
// \brief: Macro to filter tree from task output and save output trees in root files for ML studies
// \main function: FilterTree4ML
//______________________________________________________________________________________________

//______________________________________________________________________________________________
void FilterTree4ML(TString cfgFileName="config_skim_Dplus_pp5TeV.yml")
{
    const char bitSignal = 0x01;
    const char bitBkg    = 0x02;
    const char bitPrompt = 0x04;
    const char bitFD     = 0x08;
    const char bitRefl   = 0x10;

    //Load configs from yaml file
    YAML::Node config = YAML::LoadFile(cfgFileName.Data());
    if (config.IsNull()) {
        cerr << "Yaml config file not found! Exit" << endl;
        return;
    }

    const vector<string> inFileNames = config["infile"]["filename"].as<vector<string> >();
    string inDirName = config["infile"]["dirname"].as<string>();
    string inTreeName = config["infile"]["treename"].as<string>();
    bool isMC = static_cast<bool>(config["infile"]["isMC"].as<int>());

    string outTreeName = config["outfile"]["treename"].as<string>();
    string outSuffix = config["outfile"]["suffix"].as<string>();
    string outDirName = config["outfile"]["dirpath"].as<string>();

    vector<string> colsToKeep = config["skimming"]["colstokeep"].as<vector<string> >();
    if(std::find(colsToKeep.begin(), colsToKeep.end(), "inv_mass") == colsToKeep.end())
        cout << "Warning: invariant mass branch (inv_mass) disabled. Are you sure you don't want to keep it?" << endl;
    if(std::find(colsToKeep.begin(), colsToKeep.end(), "pt_cand") == colsToKeep.end())
        cout << "Warning: pt branch (pt_cand) disabled. Are you sure you don't want to keep it?" << endl;

    double PtMin = config["skimming"]["pt"]["min"].as<double>();
    double PtMax = config["skimming"]["pt"]["max"].as<double>();

    ROOT::EnableImplicitMT(); //tell ROOT to go parallel
    ROOT::RDataFrame dataFrame(Form("%s/%s", inDirName.data(), inTreeName.data()), inFileNames);
    
    if(colsToKeep.size() == 0)
    {
        colsToKeep = dataFrame.GetColumnNames();
        std::remove(colsToKeep.begin(), colsToKeep.end(), "cand_type"); //just remove cand_type if not explicitly kept
        colsToKeep.pop_back();
    }

    //select desired pT bin
    cout << "Applying pT selection" << endl;
    auto dataFramePtCut = dataFrame.Filter(Form("pt_cand > %f && pt_cand < %f", PtMin, PtMax));

    if(isMC)
    {
        cout << "Getting bkg dataframe" << endl;
        auto dataFramePtCutBkg = dataFrame.Filter(Form("cand_type >> %c & 0x01", bitBkg));
        cout << "Getting prompt dataframe" << endl;
        auto dataFramePtCutPrompt = dataFrame.Filter(Form("(cand_type >> %c & 0x01) && (cand_type >> %c & 0x01) && !(cand_type >> %c & 0x01)", bitSignal, bitPrompt, bitRefl));
        cout << "Getting FD dataframe" << endl;
        auto dataFramePtCutFD = dataFrame.Filter(Form("(cand_type >> %c & 0x01) && (cand_type >> %c & 0x01) && !(cand_type >> %c & 0x01)", bitSignal, bitFD, bitRefl));
        cout << "Getting reflected prompt dataframe" << endl;
        auto dataFramePtCutPromptRefl = dataFrame.Filter(Form("(cand_type >> %c & 0x01) && (cand_type >> %c & 0x01) && (cand_type >> %c & 0x01)", bitSignal, bitPrompt, bitRefl));
        cout << "Getting reflected signal dataframe" << endl;
        auto dataFramePtCutFDRefl = dataFrame.Filter(Form("(cand_type >> %c & 0x01) && (cand_type >> %c & 0x01) && (cand_type >> %c & 0x01)", bitSignal, bitFD, bitRefl));

        if(*dataFramePtCutBkg.Count() > 0)
        {
            cout << "Saving bkg tree" << endl;
            dataFramePtCutBkg.Snapshot(outTreeName.data(), Form("%s/Bkg%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax));
        }
        if(*dataFramePtCutPrompt.Count() > 0)
        {
            cout << "Saving prompt tree" << endl;
            dataFramePtCutPrompt.Snapshot(outTreeName.data(), Form("%s/Prompt%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax));
        }
        if(*dataFramePtCutFD.Count() > 0)
        {
            cout << "Saving FD tree" << endl;
            dataFramePtCutFD.Snapshot(outTreeName.data(), Form("%s/FD%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax));
        }
        if(*dataFramePtCutPromptRefl.Count() > 0)
        {    
            cout << "Saving prompt reflected tree" << endl;
            dataFramePtCutPromptRefl.Snapshot(outTreeName.data(), Form("%s/PromptRefl%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax));
        }
        if(*dataFramePtCutFDRefl.Count() > 0)
        {
            cout << "Saving FD reflected tree" << endl;
            dataFramePtCutFDRefl.Snapshot(outTreeName.data(), Form("%s/FDRefl%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax));
        }
    }
    else
    {
        cout << "Saving data tree" << endl;
        dataFramePtCut.Snapshot(outTreeName.data(), Form("%s/Data%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax), colsToKeep);
    }
}

    