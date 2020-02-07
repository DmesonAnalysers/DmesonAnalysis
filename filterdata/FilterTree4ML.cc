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
    const int bitSignal  = BIT(0);
    const int bitBkg     = BIT(1);
    const int bitPrompt  = BIT(2);
    const int bitFD      = BIT(3);
    const int bitRefl    = BIT(4);
    const int bitSecPeak = BIT(9);

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

    string preSelections = "";
    if(config["skimming"]["preselections"].Type() != YAML::NodeType::Null)
        preSelections = config["skimming"]["preselections"].as<string>();

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
    cout << "Applying selections" << endl;
    TString totsel = Form("pt_cand > %f && pt_cand < %f", PtMin, PtMax);
    if(preSelections != "")
    {
        totsel.Append(" & ");
        totsel.Append(preSelections.data());
    }
    auto dataFramePtCutSel = dataFrame.Filter(totsel.Data());

    if(isMC)
    {
        cout << "Getting bkg dataframe" << endl;
        auto dataFramePtCutSelBkg = dataFramePtCutSel.Filter(Form("(cand_type & %d) > 0", bitBkg));
        cout << "Getting prompt dataframe" << endl;
        auto dataFramePtCutSelPrompt = dataFramePtCutSel.Filter(Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitPrompt, bitRefl));
        cout << "Getting FD dataframe" << endl;
        auto dataFramePtCutSelFD = dataFramePtCutSel.Filter(Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitFD, bitRefl));
        cout << "Getting reflected prompt dataframe" << endl;
        auto dataFramePtCutSelPromptRefl = dataFramePtCutSel.Filter(Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0", bitSignal, bitPrompt, bitRefl));
        cout << "Getting reflected signal dataframe" << endl;
        auto dataFramePtCutSelFDRefl = dataFramePtCutSel.Filter(Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0", bitSignal, bitFD, bitRefl));
        cout << "Getting second-peak prompt dataframe" << endl;
        auto dataFramePtCutSelSecPeakPrompt = dataFramePtCutSel.Filter(Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSecPeak, bitPrompt, bitRefl));
        cout << "Getting second-peak FD dataframe" << endl;
        auto dataFramePtCutSelSecPeakFD = dataFramePtCutSel.Filter(Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSecPeak, bitFD, bitRefl));

        if(*dataFramePtCutSelBkg.Count() > 0)
        {
            cout << "Saving bkg tree" << endl;
            dataFramePtCutSelBkg.Snapshot(outTreeName.data(), Form("%s/Bkg%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax), colsToKeep);
        }
        if(*dataFramePtCutSelPrompt.Count() > 0)
        {
            cout << "Saving prompt tree" << endl;
            dataFramePtCutSelPrompt.Snapshot(outTreeName.data(), Form("%s/Prompt%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax), colsToKeep);
        }
        if(*dataFramePtCutSelFD.Count() > 0)
        {
            cout << "Saving FD tree" << endl;
            dataFramePtCutSelFD.Snapshot(outTreeName.data(), Form("%s/FD%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax), colsToKeep);
        }
        if(*dataFramePtCutSelPromptRefl.Count() > 0)
        {    
            cout << "Saving prompt reflected tree" << endl;
            dataFramePtCutSelPromptRefl.Snapshot(outTreeName.data(), Form("%s/PromptRefl%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax), colsToKeep);
        }
        if(*dataFramePtCutSelFDRefl.Count() > 0)
        {
            cout << "Saving FD reflected tree" << endl;
            dataFramePtCutSelFDRefl.Snapshot(outTreeName.data(), Form("%s/FDRefl%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax), colsToKeep);
        }
        if(*dataFramePtCutSelSecPeakPrompt.Count() > 0)
        {
            cout << "Saving prompt tree" << endl;
            dataFramePtCutSelSecPeakPrompt.Snapshot(outTreeName.data(), Form("%s/SecPeakPrompt%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax), colsToKeep);
        }
        if(*dataFramePtCutSelSecPeakFD.Count() > 0)
        {
            cout << "Saving FD tree" << endl;
            dataFramePtCutSelSecPeakFD.Snapshot(outTreeName.data(), Form("%s/SecPeakFD%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax), colsToKeep);
        }
    }
    else
    {
        cout << "Saving data tree" << endl;
        dataFramePtCutSel.Snapshot(outTreeName.data(), Form("%s/Data%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin, PtMax), colsToKeep);
    }
}

    