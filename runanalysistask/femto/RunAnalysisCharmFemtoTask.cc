#if !defined(__CINT__) || defined(__CLING__)

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
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
#include "AliAnalysisTaskCharmingFemto.h"
#include "AliAnalysisTaskSECharmHadronMLSelector.h"

#endif

using namespace std;

//______________________________________________
void RunAnalysisCharmFemtoTask(TString configfilename, TString runMode = "full", bool mergeviajdl = true)
{

    //_________________________________________________________________________________________________________________
    //load config
    YAML::Node config = YAML::LoadFile(configfilename.Data());
    string aliPhysVersion = config["AliPhysicVersion"].as<string>();
    string gridDataDir = config["datadir"].as<string>();
    string gridDataPattern = config["datapattern"].as<string>();
    string gridWorkingDir = config["gridworkdir"].as<string>();

    bool useMCTruthReco = static_cast<bool>(config["useMCTruthReco"].as<int>());
    bool useMCTruthGen = static_cast<bool>(config["useMCTruthGen"].as<int>());
    bool useTree = config["task"]["useTree"].as<bool>();
    string cutVariation = config["cutvariation"].as<string>();

    int splitmaxinputfilenum = config["splitmaxinputfilenum"].as<int>();
    int nmasterjobs = config["nmasterjobs"].as<int>();

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

    //task options
    string cutFileName = config["task"]["cuts"]["infile"].as<string>();
    string cutObjName = config["task"]["cuts"]["objname"].as<string>();
    string triggerMask = config["task"]["triggermask"].as<string>();
    bool applyML = static_cast<bool>(config["task"]["applyML"]["doapplyML"].as<int>());
    string confFileML = config["task"]["applyML"]["configfile"].as<string>();
    int pdgLight = config["task"]["pdglight"].as<int>();
    std::string charmDecChannel = config["task"]["charmDecChannel"].as<std::string>();
    std::string massSelection = config["task"]["massselection"].as<std::string>();

    bool useMLselectorTask = static_cast<bool>(config["task"]["applyML"]["MLselector"]["enable"].as<int>());
    string MLSelcutFileName = config["task"]["applyML"]["MLselector"]["infile"].as<string>();
    string MLSelcutObjName = config["task"]["applyML"]["MLselector"]["objname"].as<string>();
    string MLSelconfFileML = config["task"]["applyML"]["MLselector"]["configfile"].as<string>();

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

    //mult selection task
    AliMultSelectionTask *multSel = reinterpret_cast<AliMultSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"))));

    //pid response task
    AliAnalysisTaskPIDResponse *pidResp = reinterpret_cast<AliAnalysisTaskPIDResponse *>(gInterpreter->ProcessLine(Form(".x %s(%d)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"), isRunOnMC)));

    AliAnalysisTaskSECharmHadronMLSelector *taskMLSel = nullptr;
    if(useMLselectorTask)
    {
        std::string charmDecChannelMLSelector = charmDecChannel;
        if (charmDecChannel == "kDstartoKpipi") {
            charmDecChannelMLSelector = "kDstartoD0pi";
        }
        taskMLSel = reinterpret_cast<AliAnalysisTaskSECharmHadronMLSelector*>(gInterpreter->ProcessLine(Form(".x %s(\"%s\", \"%s\", AliAnalysisTaskSECharmHadronMLSelector::%s, \"%s\", \"""\", %s)", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskCharmHadronMLSelector.C"), MLSelcutFileName.data(), MLSelconfFileML.data(), charmDecChannelMLSelector.data(), MLSelcutObjName.data(), triggerMask.data())));
    }
    AliAnalysisTaskCharmingFemto *task = reinterpret_cast<AliAnalysisTaskCharmingFemto*>(gInterpreter->ProcessLine(Form(".x %s(%d, %d, %d, %d, true, \"%s\", AliAnalysisTaskCharmingFemto::%s, \"%s\", \"%s\", \"""\", %d, \"%s\", 0, AliAnalysisTaskCharmingFemto::%s, %d, \"%s\")", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/macros/AddTaskAnyCharmingFemto.C"), isRunOnMC, useMCTruthReco, useMCTruthGen, useTree, triggerMask.data(), charmDecChannel.data(), cutFileName.data(), cutObjName.data(), applyML, confFileML.data(), massSelection.data(), pdgLight, cutVariation.data())));
    if(useMLselectorTask)
        task->SetIsDependentOnMLSelector();

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
        for (const auto &path : std::filesystem::recursive_directory_iterator(pathToLocalAODfiles.data())) {
            auto currentPath = path.path();
            if (std::string(currentPath).find("AliAOD.root") != std::string::npos) {
                std::cout << "Adding data file:      " << currentPath << '\n';
                chainAOD->Add(std::string(currentPath).data());
            } else if (std::string(currentPath).find("AliAOD.VertexingHF.root") != std::string::npos) {
                std::cout << "Adding vertexing file: " << currentPath << '\n';
                chainAODfriend->Add(std::string(currentPath).data());
            }
        }

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
        alienHandler->SetNrunsPerMaster(nruns/nmasterjobs);

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
