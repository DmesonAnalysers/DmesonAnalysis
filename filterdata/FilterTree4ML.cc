
#if !defined (__CINT__) || defined (__CLING__)

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include "yaml-cpp/yaml.h"

#include <TFile.h>
#include <TTree.h>
#include "ROOT/RDataFrame.hxx"

#endif

using std::vector;
using std::map;
using std::string;
using std::cout;
using std::cerr;
using std::endl;

//______________________________________________________________________________________________
// \brief: Macro to filter tree from task output and save output trees in root files for ML studies
// \main function: FilterTree4ML
//______________________________________________________________________________________________

//______________________________________________________________________________________________
void FilterTree4ML(TString cfgFileName="config_Dstar_data_skim_pp5TeV.yml")
{
    // Common bits
    const int bitSignal       = BIT(0);
    const int bitBkg          = BIT(1);
    const int bitPrompt       = BIT(2);
    const int bitFD           = BIT(3);
    const int bitRefl         = BIT(4);
    // Channel specific bits
    // Ds
    const int bitSecPeakDs    = BIT(9);
    // LctopK0s
    const int bitLctopK0s     = BIT(9);
    // LctopiL
    const int bitLctopiL      = BIT(10);
    // LctopKpi
    const int bitLcNonRes     = BIT(9);
    const int bitLcLambda1520 = BIT(10);
    const int bitLcKStar      = BIT(11);
    const int bitLcDelta      = BIT(12);

    // Load configs from yaml file
    YAML::Node config = YAML::LoadFile(cfgFileName.Data());
    if (config.IsNull()) {
        cerr << "Error: Yaml config file not found! Exit" << endl;
        return;
    }

    vector<string> channels = {"Ds", "D0", "Dplus", "Dstar", "LctopKpi", "LctopK0s", "LctopLi"};
    string channel = config["channel"].as<string>();
    if (std::find(channels.begin(), channels.end(), channel) == channels.end()) {
        cerr << "Error: only Ds, D0, Dplus, LctopKpi, LctopK0s, and LctopiL channels are implemented! Exit" << endl;
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

    const vector<double> PtMin = config["skimming"]["pt"]["min"].as<vector<double>>();
    const vector<double> PtMax = config["skimming"]["pt"]["max"].as<vector<double>>();

    ROOT::EnableImplicitMT(32); // tell ROOT to go parallel

    for(uint iPtBin = 0; iPtBin < PtMin.size(); iPtBin++){
        ROOT::RDataFrame dataFrame(Form("%s/%s", inDirName.data(), inTreeName.data()), inFileNames);

        if(colsToKeep.size() == 0)
        {
            colsToKeep = dataFrame.GetColumnNames();
            std::remove(colsToKeep.begin(), colsToKeep.end(), "cand_type"); // just remove cand_type if not explicitly kept
            colsToKeep.pop_back();
        }

        // Select desired pT bin
        cout << "Applying selections" << endl;
        TString totsel = Form("pt_cand > %f && pt_cand < %f", PtMin[iPtBin], PtMax[iPtBin]);
        if(preSelections != "")
        {
            totsel.Append(" & ");
            totsel.Append(preSelections.data());
        }
        auto dataFramePtCutSel = dataFrame.Filter(totsel.Data());

        // Add single track variables as those used in AOD filtering
        auto ptMinFormula = [](float pt0, float pt1, float pt2)
        {
            vector<float> ptVec = {pt0, pt1, pt2};
            return *std::min_element(ptVec.begin(), ptVec.end());
        };
        auto d0MinFormula = [](float pt0, float pt1, float pt2, float d00, float d01, float d02)
        {
            vector<float> ptVec = {pt0, pt1, pt2};
            vector<float> d0Vec = {fabs(d00), fabs(d01), fabs(d02)};
            vector<float> d0VecSel{};
            for(unsigned int iDau=0; iDau<ptVec.size(); iDau++)
            {
                if(ptVec[iDau] < 2)
                    d0VecSel.push_back(d0Vec[iDau]);
            }
            if(d0VecSel.size() == 0)
                return 999.f;

            return *std::min_element(d0Vec.begin(), d0Vec.end());
        };

        if(config["singletrackvars"]["addAODfiltervars"] && config["singletrackvars"]["addAODfiltervars"].as<int>())
        {
            colsToKeep.push_back("pt_prong_min");
            colsToKeep.push_back("imp_par_min_ptgtr2");
        }

        std::map<string, string> labelsContr = {{"bkg", "Bkg"}, {"prompt_sig", "Prompt"}, {"FD_sig", "FD"},
                                                {"prompt_sig_refl", "PromptRefl"}, {"FD_sig_refl", "FDRefl"},
                                                {"prompt_sec_peak", "SecPeakPrompt"}, {"FD_sec_peak", "SecPeakFD"},
                                                {"prompt_sig_nonreso", "PromptNonRes"}, {"FD_sig_nonreso", "FDNonRes"},
                                                {"prompt_sig_Lambda1520", "PromptLambda1520"}, {"FD_sig_Lambda1520", "FDLambda1520"},
                                                {"prompt_sig_KStar", "PromptKStar"}, {"FD_sig_KStar", "FDKStar"},
                                                {"prompt_sig_Delta", "PromptDelta"}, {"FD_sig_Delta", "FDDelta"}};
        std::map<string, string> bitsForSel{};

        if(isMC)
        {
            if(channel == "Ds") {
                bitsForSel["bkg"]                   = Form("(cand_type & %d) > 0", bitBkg);
                bitsForSel["prompt_sig"]            = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitPrompt, bitRefl);
                bitsForSel["FD_sig"]                = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitFD, bitRefl);
                bitsForSel["prompt_sig_refl"]       = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0", bitSignal, bitPrompt, bitRefl);
                bitsForSel["FD_sig_refl"]           = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0", bitSignal, bitFD, bitRefl);
                bitsForSel["prompt_sec_peak"]       = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSecPeakDs, bitPrompt, bitRefl);
                bitsForSel["FD_sec_peak"]           = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSecPeakDs, bitFD, bitRefl);
            }
            else if(channel == "Dplus") {
                bitsForSel["bkg"]                   = Form("(cand_type & %d) > 0", bitBkg);
                bitsForSel["prompt_sig"]            = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitPrompt, bitRefl);
                bitsForSel["FD_sig"]                = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitFD, bitRefl);
            }
            else if(channel == "D0") {
                bitsForSel["bkg"]                   = Form("(cand_type & %d) > 0", bitBkg);
                bitsForSel["prompt_sig"]            = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitPrompt, bitRefl);
                bitsForSel["FD_sig"]                = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitFD, bitRefl);
                bitsForSel["prompt_sig_refl"]       = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0", bitSignal, bitPrompt, bitRefl);
                bitsForSel["FD_sig_refl"]           = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0", bitSignal, bitFD, bitRefl);
            }
             else if(channel == "Dstar") {
                bitsForSel["bkg"]                   = Form("(cand_type & %d) > 0", bitBkg);
                bitsForSel["prompt_sig"]            = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitPrompt, bitRefl);
                bitsForSel["FD_sig"]                = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitFD, bitRefl);
            }
            else if(channel == "LctopKpi") {
                bitsForSel["bkg"]                   = Form("(cand_type & %d) > 0", bitBkg);
                bitsForSel["prompt_sig_nonreso"]    = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLcNonRes, bitPrompt, bitRefl);
                bitsForSel["FD_sig_nonreso"]        = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLcNonRes, bitFD, bitRefl);
                bitsForSel["prompt_sig_Lambda1520"] = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLcLambda1520, bitPrompt, bitRefl);
                bitsForSel["FD_sig_Lambda1520"]     = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLcLambda1520, bitFD, bitRefl);
                bitsForSel["prompt_sig_KStar"]      = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLcKStar, bitPrompt, bitRefl);
                bitsForSel["FD_sig_KStar"]          = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLcKStar, bitFD, bitRefl);
                bitsForSel["prompt_sig_Delta"]      = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLcDelta, bitPrompt, bitRefl);
                bitsForSel["FD_sig_Delta"]          = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLcDelta, bitFD, bitRefl);
                bitsForSel["prompt_sig_refl"]       = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0", bitSignal, bitPrompt, bitRefl);
                bitsForSel["FD_sig_refl"]           = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0", bitSignal, bitFD, bitRefl);
            }
            else if(channel == "LctopK0s") {
                bitsForSel["bkg"]                   = Form("(cand_type & %d) > 0", bitBkg);
                bitsForSel["prompt_sig"]            = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLctopK0s, bitPrompt, bitRefl);
                bitsForSel["FD_sig"]                = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLctopK0s, bitFD, bitRefl);
            }
            else if(channel == "LctopiL") {
                bitsForSel["bkg"]                   = Form("(cand_type & %d) > 0", bitBkg);
                bitsForSel["prompt_sig"]            = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLctopiL, bitPrompt, bitRefl);
                bitsForSel["FD_sig"]                = Form("(cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) > 0 && (cand_type & %d) == 0", bitSignal, bitLctopiL, bitFD, bitRefl);
            }

            for(auto &contr: bitsForSel) {
                cout << Form("Getting %s dataframe", labelsContr[contr.first].data()) << endl;
                auto dataFramePtCutSelContr = dataFramePtCutSel.Filter(contr.second.data());
                if(*dataFramePtCutSelContr.Count() > 0)
                {
                    cout << "Saving bkg tree" << endl;
                    if(config["singletrackvars"]["addAODfiltervars"] && config["singletrackvars"]["addAODfiltervars"].as<int>())
                        dataFramePtCutSelContr.Define("pt_prong_min", ptMinFormula, {"pt_prong0", "pt_prong1", "pt_prong2"})
                                            .Define("imp_par_min_ptgtr2", d0MinFormula, {"pt_prong0", "pt_prong1", "pt_prong2", "imp_par_prong0", "imp_par_prong1", "imp_par_prong2"})
                                            .Snapshot(outTreeName.data(), Form("%s/%s%s_pT_%0.f_%0.f.root", outDirName.data(), labelsContr[contr.first].data(), outSuffix.data(), PtMin[iPtBin], PtMax[iPtBin]), colsToKeep);
                    else
                        dataFramePtCutSelContr.Snapshot(outTreeName.data(), Form("%s/%s%s_pT_%0.f_%0.f.root", outDirName.data(), labelsContr[contr.first].data(), outSuffix.data(), PtMin[iPtBin], PtMax[iPtBin]), colsToKeep);
                }
            }
        }
        else
        {
            cout << "Saving data tree" << endl;
            if(config["singletrackvars"]["addAODfiltervars"] && config["singletrackvars"]["addAODfiltervars"].as<int>())
                dataFramePtCutSel.Define("pt_prong_min", ptMinFormula, {"pt_prong0", "pt_prong1", "pt_prong2"})
                                .Define("imp_par_min_ptgtr2", d0MinFormula, {"pt_prong0", "pt_prong1", "pt_prong2", "imp_par_prong0", "imp_par_prong1", "imp_par_prong2"})
                                .Snapshot(outTreeName.data(), Form("%s/Data%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin[iPtBin], PtMax[iPtBin]), colsToKeep);
            else
                dataFramePtCutSel.Snapshot(outTreeName.data(), Form("%s/Data%s_pT_%0.f_%0.f.root", outDirName.data(), outSuffix.data(), PtMin[iPtBin], PtMax[iPtBin]), colsToKeep);
        }
    }
}
