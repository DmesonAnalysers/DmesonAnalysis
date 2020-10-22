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
#include "AliAODInputHandler.h"

#include "AliPhysicsSelectionTask.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliMultSelectionTask.h"
#include "AliAnalysisTaskSEImproveITSCVMFS.h"
#include "AliAnalysisTaskSEImproveITS3.h"
#include "AliAnalysisTaskSECleanupVertexingHF.h"
#include "AliAnalysisTaskSEDplus.h"
#include "AliAnalysisTaskSEDs.h"
#include "AliAnalysisTaskSENonPromptLc.h"

#endif

using namespace std;

enum improverType {
    kNoImprover, 
    kCurrentImprover,
    kUpgradeImprover
};

//______________________________________________
void RunAnalysisDplusDsLcTask(TString configfilename, TString runMode = "full", bool mergeviajdl = true)
{

    //_________________________________________________________________________________________________________________
    //load config
    YAML::Node config = YAML::LoadFile(configfilename.Data());
    string aliPhysVersion = config["AliPhysicVersion"].as<string>();
    string gridDataDir = config["datadir"].as<string>();
    string gridDataPattern = config["datapattern"].as<string>();
    string gridWorkingDir = config["gridworkdir"].as<string>();
    int splitmaxinputfilenum = config["splitmaxinputfilenum"].as<int>();

    string hadron = config["hadron"].as<string>();
    string sSystem = config["system"].as<string>();
    int system = -1;
    if (sSystem == "pp")
    {
        if (hadron == "Dplus")
            system = AliAnalysisTaskSEDplus::kpp;
        else if (hadron == "Ds")
            system = AliAnalysisTaskSEDs::kpp;
    }
    else if (sSystem == "PbPb")
    {
        if (hadron == "Dplus")
            system = AliAnalysisTaskSEDplus::kPbPb;
        else if (hadron == "Ds")
            system = AliAnalysisTaskSEDs::kPbPb;
    }
    else if (sSystem == "Upgr")
    {
        if (hadron == "Ds")
            system = AliAnalysisTaskSEDs::kUpgr;
    }
    else
    {
        cerr << "ERROR: Only pp, PbPb, and Upgr systems are supported, please check your config file. Exit" << endl;
        return;
    }

    vector<int> runlist = config["runs"].as<vector<int>>();
    const int nruns = runlist.size();
    bool isRunOnMC = strstr(gridDataDir.c_str(), "sim");

    bool local = false;
    bool gridTest = false;
    string pathToLocalAODfiles = "";
    if (!gGrid)
        TGrid::Connect("alien://");
    if (config["runtype"].as<string>() == "local")
    {
        local = true;
        pathToLocalAODfiles = config["pathtolocalAOD"].as<string>();
    }
    else
    {
        if (config["runtype"].as<string>() == "test")
            gridTest = true;
    }

    //improver options
    int useImprover = kNoImprover;
    string improverPeriod = "";
    string improverFileCurrent = "";
    string improverFileUpgrade = "";
    if(config["improver"]["current"]["enable"].as<int>())
    {
        useImprover = kCurrentImprover;
        improverPeriod = config["improver"]["current"]["period"].as<string>();
    }
    if(config["improver"]["upgrade"]["enable"].as<int>())
    {
        if(useImprover == kCurrentImprover)
        {
            cerr << "ERROR: you cannot use the improver for the current ITS and for the upgraded one! Check your config file. Exit" << endl;
            return;
        }
        else
        {
            useImprover = kUpgradeImprover;
            improverFileCurrent = config["improver"]["upgrade"]["currentresolfile"].as<string>();
            improverFileUpgrade = config["improver"]["upgrade"]["upgraderesolfile"].as<string>();
        }
    }


    //task options
    string wagonName = config["task"]["wagonname"].as<string>();
    string cutFileName = config["task"]["cuts"]["infile"].as<string>();
    string cutObjName = config["task"]["cuts"]["objname"].as<string>();
    bool storeSparse = static_cast<bool>(config["task"]["storesparse"].as<int>());
    bool storeTreeML = static_cast<bool>(config["task"]["treeML"]["storetree"].as<int>());
    bool fillOnlySignalTreeML = static_cast<bool>(config["task"]["treeML"]["fillonlysignal"].as<int>());
    bool enableTrackVarsTreeML = static_cast<bool>(config["task"]["treeML"]["enabletrackvars"].as<int>());
    string sPidTreeOpt = config["task"]["treeML"]["PIDoption"].as<string>();
    int pidTreeOpt = -1;
    if(sPidTreeOpt == "kNoPID")
        pidTreeOpt = AliHFMLVarHandler::kNoPID;
    else if(sPidTreeOpt == "kNoPID")
        pidTreeOpt = AliHFMLVarHandler::kRawPID;
    else if(sPidTreeOpt == "kNsigmaPID")
        pidTreeOpt = AliHFMLVarHandler::kNsigmaPID;
    else if(sPidTreeOpt == "kNsigmaCombPID")
        pidTreeOpt = AliHFMLVarHandler::kNsigmaCombPID;
    else if(sPidTreeOpt == "kNsigmaDetAndCombPID")
        pidTreeOpt = AliHFMLVarHandler::kNsigmaDetAndCombPID;
    else if(sPidTreeOpt == "kRawAndNsigmaPID")
        pidTreeOpt = AliHFMLVarHandler::kRawAndNsigmaPID;
    bool applyML = static_cast<bool>(config["task"]["applyML"]["doapplyML"].as<int>());
    string confFileML = config["task"]["applyML"]["configfile"].as<string>();
    int nMLbins = config["task"]["applyML"]["nbins"].as<int>();
    double MLmin = config["task"]["applyML"]["min"].as<double>();
    double MLmax = config["task"]["applyML"]["max"].as<double>();
    //_________________________________________________________________________________________________________________

    // if compile a class, tell root where to look for headers
    // gInterpreter->ProcessLine(".include $ROOTSYS/include");
    // gInterpreter->ProcessLine(".include $ALICE_ROOT/include");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    //physics selection task
    AliPhysicsSelectionTask *physseltask = nullptr;
    if (runMode != "terminate")
        physseltask = reinterpret_cast<AliPhysicsSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s(%d, %d)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"), isRunOnMC, true)));

    //pid response task
    AliAnalysisTaskPIDResponse *pidResp = reinterpret_cast<AliAnalysisTaskPIDResponse *>(gInterpreter->ProcessLine(Form(".x %s(%d)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"), isRunOnMC)));

    //mult sel task
    if(system == AliAnalysisTaskSEDplus::kPbPb || system == AliAnalysisTaskSEDs::kPbPb) {
        AliMultSelectionTask *multSel = reinterpret_cast<AliMultSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"))));
        if(strstr(gridDataDir.c_str(),"LHC15o") || strstr(gridDataDir.c_str(),"LHC16i2")) multSel->SetAlternateOADBforEstimators("LHC15o-DefaultMC-HIJING");
        else if((strstr(gridDataDir.c_str(),"LHC18q") || strstr(gridDataDir.c_str(),"LHC18r") || strstr(gridDataDir.c_str(),"LHC19c")) && isRunOnMC) multSel->SetAlternateOADBFullManualBypassMC("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/data/OADB-LHC18q-DefaultMC-HIJING.root");
    }

    //improver task (if MC)
    if (isRunOnMC)
    {
        if(useImprover==kCurrentImprover && improverPeriod != "")
        {    
            AliAnalysisTaskSEImproveITSCVMFS *taskimpr = reinterpret_cast<AliAnalysisTaskSEImproveITSCVMFS *>(gInterpreter->ProcessLine(Form(".x %s(%d,\"%s\",\"%s\",%d)", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskImproveITSCVMFS.C"), false,"","",0)));
        }
        else if(useImprover==kUpgradeImprover && improverFileCurrent != "" && improverFileUpgrade != "")
        {
            AliAnalysisTaskSEImproveITS3* taskimpr = reinterpret_cast<AliAnalysisTaskSEImproveITS3 *>(gInterpreter->ProcessLine(Form(".x %s(%d,\"%s\",\"%s\",%d)", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/upgrade/AddTaskImproverUpgrade.C"), false, improverFileCurrent.data() ,improverFileUpgrade.data(), 0)));
        }
    }

    //D+ or Ds or Lc tasks
    if (hadron == "Dplus")
    {
        AliAnalysisTaskSEDplus *taskDplus = reinterpret_cast<AliAnalysisTaskSEDplus *>(gInterpreter->ProcessLine(Form(".x %s(%d,%f,%f,%d,%d,%d,%d,\"%s\",\"%s\",\"%s\")", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskDplus.C"), system, 0., 100., storeTreeML, storeSparse, false, isRunOnMC, wagonName.data(), cutFileName.data(), cutObjName.data())));
        if(applyML)
        {
            taskDplus->SetDoMLApplication();
            taskDplus->SetMLConfigFile(confFileML.data());   
            taskDplus->SetMLBinsForSparse(nMLbins, MLmin, MLmax);     
        }
        if(storeTreeML)
        {
            taskDplus->SetMLTreePIDopt(pidTreeOpt);
            taskDplus->SetMLTreeAddTrackVar(enableTrackVarsTreeML);
            taskDplus->SetFillOnlySignalInMLtree(fillOnlySignalTreeML);
        }
    }
    else if (hadron == "Ds")
    {
        AliAnalysisTaskSEDs *taskDs = reinterpret_cast<AliAnalysisTaskSEDs *>(gInterpreter->ProcessLine(Form(".x %s(%d,%d,%d,%d,\"%s\",\"%s\",%d,%d,\"%s\",\"%s\",%d)", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskDs.C"), system, isRunOnMC, 0, storeSparse, cutFileName.data(), wagonName.data(), storeTreeML, applyML, confFileML.data(), cutObjName.data(), storeSparse)));
        if(storeSparse)
        {
            taskDs->SetFillNSparseMCAccWoQuark();
        }
        if(applyML)
        {
            taskDs->SetDoMLApplication();
            taskDs->SetMLConfigFile(confFileML.data());   
            taskDs->SetMLBinsForSparse(nMLbins, MLmin, MLmax);     
        }
        if(storeTreeML)
        {
            taskDs->SetMLTreePIDopt(pidTreeOpt);
            taskDs->SetMLTreeAddTrackVar(enableTrackVarsTreeML);
            taskDs->SetFillOnlySignalInMLtree(fillOnlySignalTreeML);
        }
    }
    else if (hadron == "LcpKpi" || hadron == "LcpK0s" || hadron == "LcpiL")
    {
        int channel = AliAnalysisTaskSENonPromptLc::kLctopKpi;
        if(hadron == "LcpK0s")
            channel = AliAnalysisTaskSENonPromptLc::kLctopK0s;
        else if(hadron == "LcpiL")
            channel = AliAnalysisTaskSENonPromptLc::kLctopiL;
        
        AliAnalysisTaskSENonPromptLc *taskLc = reinterpret_cast<AliAnalysisTaskSENonPromptLc *>(gInterpreter->ProcessLine(Form(".x %s(%d, %d, \"%s\", \"%s\", %d, \"%s\")", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskNonPromptLc.C"), channel, isRunOnMC, cutFileName.data(), wagonName.data(), storeTreeML, "AnalysisCuts")));
        if(storeTreeML)
        {
            taskLc->SetMLTreePIDopt(pidTreeOpt);
            taskLc->SetMLTreeAddTrackVar(enableTrackVarsTreeML);
            taskLc->SetFillOnlySignalInMLtree(fillOnlySignalTreeML);
        }
    }

    //cleanup task (if PbPb)
    if (system == AliAnalysisTaskSEDplus::kPbPb || system == AliAnalysisTaskSEDs::kPbPb || system == AliAnalysisTaskSEDs::kUpgr)
    {
        AliAnalysisTaskSECleanupVertexingHF *taskclean = reinterpret_cast<AliAnalysisTaskSECleanupVertexingHF *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskCleanupVertexingHF.C"))));
    }

    if (!mgr->InitAnalysis())
        return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if (local)
    {
        // if you want to run locally, we need to define some input
        TChain *chainAOD = new TChain("aodTree");
        TChain *chainAODfriend = new TChain("aodTree");

        // add a few files to the chain (change this so that your local files are added)
        chainAOD->Add(Form("%s/AliAOD.root", pathToLocalAODfiles.data()));
        chainAODfriend->Add(Form("%s/AliAOD.VertexingHF.root", pathToLocalAODfiles.data()));

        chainAOD->AddFriend(chainAODfriend);

        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chainAOD);
    }
    else
    {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();

        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

        //make sure your source files get copied to grid
        // alienHandler->SetAdditionalLibs("AliAnalysisTaskSEDplus.h AliAnalysisTaskSEDplus.cxx");
        // alienHandler->SetAnalysisSource("AliAnalysisTaskSEDplus.cxx");

        // select the aliphysics version.
        alienHandler->SetAliPhysicsVersion(aliPhysVersion.data());

        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");

        // select the input data
        alienHandler->SetGridDataDir(gridDataDir.data());
        alienHandler->SetDataPattern(Form("%s/*AliAOD.root", gridDataPattern.data()));
        alienHandler->SetFriendChainName("AliAOD.VertexingHF.root");

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
        alienHandler->SetTTL(30000);
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
