#if !defined(__CINT__) || defined(__CLING__)

#include <iostream>
#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include <TChain.h>
#include <TGrid.h>
#include <TNtuple.h>
#include <TInterpreter.h>
#include <TSystem.h>

#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisGrid.h"

#include "AliAnalysisTaskCheckHFMCProd.h"

#endif

using namespace std;

//______________________________________________
void RunAnalysisCheckHFProd(TString configfilename, TString runMode = "full", bool mergeviajdl = true)
{

    //_________________________________________________________________________________________________________________
    //load config
    YAML::Node config = YAML::LoadFile(configfilename.Data());
    string aliPhysVersion = config["AliPhysicVersion"].as<string>();
    string gridDataDir = config["datadir"].as<string>();
    string gridDataPattern = config["datapattern"].as<string>();
    string gridWorkingDir = config["gridworkdir"].as<string>();
    int splitmaxinputfilenum = config["splitmaxinputfilenum"].as<int>();

    string sSystem = config["system"].as<string>();
    int system = -1;
    if (sSystem == "pp")
        system = 0;
    else if (sSystem == "pPb")
        system = 1;
    else if (sSystem == "PbPb")
        system = 2;
    else
    {
        cerr << "ERROR: Only pp, pPb, and PbPb systems are supported, please check your config file. Exit" << endl;
        return;
    }

    vector<int> runlist = config["runs"].as<vector<int>>();
    const int nruns = runlist.size();
    bool isRunOnMC = strstr(gridDataDir.c_str(), "sim");

    bool local = false;
    bool gridTest = false;
    string pathToLocalESDfiles = "";
    if (!gGrid)
        TGrid::Connect("alien://");
    if (config["runtype"].as<string>() == "local")
    {
        local = true;
        pathToLocalESDfiles = config["pathtolocalESD"].as<string>();
    }
    else
    {
        if (config["runtype"].as<string>() == "test")
            gridTest = true;
    }

    //task options
    int nPtBins = config["task"]["ptbins"]["nbins"].as<int>();
    double ptMin = config["task"]["ptbins"]["min"].as<double>();
    double ptMax = config["task"]["ptbins"]["max"].as<double>();
    //_________________________________________________________________________________________________________________

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager();
    AliVEventHandler* esdH = new AliESDInputHandler;
    AliMCEventHandler* handlerMC = new AliMCEventHandler;
    mgr->SetInputEventHandler(esdH);
    esdH->SetNeedField();
    mgr->SetMCtruthEventHandler(handlerMC);

    // CheckHFProd task
    AliAnalysisTaskCheckHFMCProd* taskHFMCProd = reinterpret_cast<AliAnalysisTaskCheckHFMCProd *>(gInterpreter->ProcessLine(Form(".x %s(%d, %d)", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddHFMCCheck.C"), system, isRunOnMC)));
    cout << ptMin << "  " << ptMax << "  " << nPtBins << endl;
    taskHFMCProd->SetPtBins(ptMin, ptMax, nPtBins);
    taskHFMCProd->SetYBins(-2., 2., 40);

    if (!mgr->InitAnalysis())
        return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if (local)
    {
        // if you want to run locally, we need to define some input
        TChain *chainESD = new TChain("esdTree");
        TChain *chainESDfriend = new TChain("esdFriendTree");

        // add a few files to the chain (change this so that your local files are added)
        chainESD->Add(Form("%s/AliESDs.root", pathToLocalESDfiles.data()));
        chainESDfriend->Add(Form("%s/AliAOD.VertexingHF.root", pathToLocalESDfiles.data()));
        chainESD->AddFriend(chainESDfriend);

        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chainESD);
    }
    else
    {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();

        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

        // select the aliphysics version.
        alienHandler->SetAliPhysicsVersion(aliPhysVersion.data());

        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");

        // select the input data
        alienHandler->SetGridDataDir(gridDataDir.data());
        alienHandler->SetDataPattern(Form("%s/*AliESDs.root", gridDataPattern.data()));

        // MC has no prefix, data has prefix 000
        if (!isRunOnMC)
            alienHandler->SetRunPrefix("000");

        for (int iRun = 0; iRun < nruns; iRun++)
            alienHandler->AddRunNumber(runlist[iRun]);
        alienHandler->SetNrunsPerMaster(1);

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(splitmaxinputfilenum);
        alienHandler->SetExecutable("myAnalysis.sh");

        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetJDLName(Form("%s.jdl", gridWorkingDir.data()));

        alienHandler->SetOutputToRunNo(true);
        alienHandler->SetKeepLogs(true);

        // merging: run with true to merge on grid
        // after re-running the jobs in SetRunMode("terminate")
        // (see below) mode, set SetMergeViaJDL(false)
        // to collect final results
        alienHandler->SetMaxMergeStages(3); //2, 3
        alienHandler->SetMergeViaJDL(mergeviajdl);

        // define the output folders
        alienHandler->SetGridWorkingDir(gridWorkingDir.data());
        alienHandler->SetGridOutputDir("output");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);

        if (gridTest)
        {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        }
        else
        {
            // else launch the full grid analysis
            alienHandler->SetRunMode(runMode.Data()); //terminate
            mgr->StartAnalysis("grid");
        }
    }
}