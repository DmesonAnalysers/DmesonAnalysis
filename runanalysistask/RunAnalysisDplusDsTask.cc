#if !defined(__CINT__) || defined(__CLING__)

#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include <TChain.h>
#include <TGrid.h>
#include <TNtuple.h>

#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"

#include "AliPhysicsSelectionTask.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliMultSelectionTask.h"
#include "AliAnalysisTaskSEImproveITS.h"
#include "AliAnalysisTaskSECleanupVertexingHF.h"
#include "AliAnalysisTaskSEDplus.h"
#include "AliAnalysisTaskSEDs.h"

#endif

using namespace std;

//______________________________________________
void RunAnalysisDplusDsTask(TString configfilename, TString runMode = "full", bool mergeviajdl = true)
{

    //_________________________________________________________________________________________________________________
    //load config
    YAML::Node config = YAML::LoadFile(configfilename.Data());
    string aliPhysVersion = config["AliPhysicVersion"].as<string>();
    string gridDataDir = config["datadir"].as<string>();
    string gridDataPattern = config["datapattern"].as<string>();
    string gridWorkingDir = config["gridworkdir"].as<string>();
    int splitmaxinputfilenum = config["splitmaxinputfilenum"].as<int>();

    string meson = config["meson"].as<string>();
    string sSystem = config["system"].as<string>();
    int system = -1;
    if (sSystem == "pp")
    {
        if (meson == "Dplus")
            system = AliAnalysisTaskSEDplus::kpp;
        else if (meson == "Ds")
            system = AliAnalysisTaskSEDs::kpp;
    }
    else if (sSystem == "PbPb")
    {
        if (meson == "Dplus")
            system = AliAnalysisTaskSEDplus::kPbPb;
        else if (meson == "Ds")
            system = AliAnalysisTaskSEDs::kPbPb;
    }
    else
    {
        cerr << "ERROR: Only pp and PbPb are supported, please check your config file. Exit" << endl;
        return;
    }

    vector<int> runlist = config["runs"].as<vector<int>>();
    const int nruns = runlist.size();
    bool isRunOnMC = strstr(gridDataDir.c_str(), "sim");

    bool local = false;
    bool gridTest = false;
    string pathToLocalAODfiles = "";
    if (config["runtype"].as<string>() == "local")
    {
        local = true;
        pathToLocalAODfiles = config["pathtolocalAOD"].as<string>();
    }
    else
    {
        if (!gGrid)
            TGrid::Connect("alien://");

        if (config["runtype"].as<string>() == "test")
            gridTest = true;
    }

    bool useImprover = static_cast<bool>(config["improver"]["enable"].as<int>());
    string improverPeriod = config["improver"]["period"].as<string>();

    string wagonName = config["wagonname"].as<string>();
    string cutFileName = config["cuts"]["infile"].as<string>();
    string cutObjName = config["cuts"]["objname"].as<string>();

    //task options
    bool storeSparse = static_cast<bool>(config["task"]["storesparse"].as<int>());
    bool storeTreeML = static_cast<bool>(config["task"]["ML"]["tree"]["storetree"].as<int>());
    bool fillOnlySignalTreeML = static_cast<bool>(config["task"]["ML"]["tree"]["fillonlysignal"].as<int>());
    bool enableTrackVarsTreeML = static_cast<bool>(config["task"]["ML"]["tree"]["enabletrackvars"].as<int>());
    string sPidTreeOpt = config["task"]["ML"]["tree"]["PIDoption"].as<string>();
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
    bool applyML = static_cast<bool>(config["task"]["ML"]["application"]["doapplyML"].as<int>());
    string confFileML = config["task"]["ML"]["application"]["configfile"].as<string>();
    int nMLbins = config["task"]["ML"]["application"]["bins4sparse"]["nbins"].as<int>();
    double MLmin = config["task"]["ML"]["application"]["bins4sparse"]["min"].as<double>();
    double MLmax = config["task"]["ML"]["application"]["bins4sparse"]["max"].as<double>();
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
        physseltask = reinterpret_cast<AliPhysicsSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s (%d)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"), isRunOnMC)));

    //pid response task
    AliAnalysisTaskPIDResponse *pidResp = reinterpret_cast<AliAnalysisTaskPIDResponse *>(gInterpreter->ProcessLine(Form(".x %s (%d)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"), isRunOnMC)));

    //mult sel task
    if(system == AliAnalysisTaskSEDplus::kPbPb || system == AliAnalysisTaskSEDs::kPbPb) {
        AliMultSelectionTask *multSel = reinterpret_cast<AliMultSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"))));

        if(strstr(gridDataDir.c_str(),"LHC15o") || strstr(gridDataDir.c_str(),"LHC16i2")) multSel->SetAlternateOADBforEstimators("LHC15o-DefaultMC-HIJING");
        else if((strstr(gridDataDir.c_str(),"LHC18q") || strstr(gridDataDir.c_str(),"LHC18r") || strstr(gridDataDir.c_str(),"LHC19c")) && isRunOnMC) multSel->SetAlternateOADBFullManualBypassMC("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/data/OADB-LHC18q-DefaultMC-HIJING.root");
    }

    //improver task (if MC)
    if (isRunOnMC && useImprover && improverPeriod != "")
    {
        AliAnalysisTaskSEImproveITS *taskimpr = reinterpret_cast<AliAnalysisTaskSEImproveITS *>(gInterpreter->ProcessLine(Form(".x %s(%d,\"%s\")", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskImprover.C"), false, improverPeriod.data())));
    }

    //D+ or Ds tasks
    if (meson == "Dplus")
    {
        AliAnalysisTaskSEDplus *taskDplus = reinterpret_cast<AliAnalysisTaskSEDplus *>(gInterpreter->ProcessLine(Form(".x %s(%d,%f,%f,%d,%d,%d,%d,\"%s\",\"%s\",\"%s\")", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/AddTaskDplus.C"), system, 0., 100., storeTreeML, storeSparse, false, isRunOnMC, wagonName.data(), cutFileName.data(), cutObjName.data())));
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
    else if (meson == "Ds")
    {
        AliAnalysisTaskSEDs *taskDs = reinterpret_cast<AliAnalysisTaskSEDs *>(gInterpreter->ProcessLine(Form(".x %s(%d,%d,%d,%d,\"%s\",\"%s\",%d,%d,\"%s\",\"%s\",%d)", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/AddTaskDs.C"), system, isRunOnMC, 1, storeSparse, cutFileName.data(), wagonName.data(), storeTreeML, applyML, confFileML.data(), cutObjName.data(), storeSparse)));
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

    //cleanup task (if PbPb)
    if (system == AliAnalysisTaskSEDplus::kPbPb || system == AliAnalysisTaskSEDs::kPbPb)
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
        alienHandler->SetExecutable(Form("%s.sh", gridDataDir.data()));

        // specify how many seconds your job may take
        alienHandler->SetTTL(30000);
        alienHandler->SetJDLName(Form("%s.jdl", gridDataDir.data()));

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
