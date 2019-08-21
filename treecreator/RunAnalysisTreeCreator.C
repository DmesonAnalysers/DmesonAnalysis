#if !defined (__CINT__) || defined (__CLING__)

#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include <TChain.h>
#include <TGrid.h>

#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"

#include "AliPhysicsSelectionTask.h"
#include "AliMultSelectionTask.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAnalysisTaskSECleanupVertexingHF.h"

#include "AliAnalysisTaskSEHFTreeCreator.h"
#include "AliAnalysisTaskSEDs.h"

#endif

using namespace std;
enum system {kpp, kpPb, kPbPb};

//______________________________________________
//SETTINGS

bool runLocal=true;                                      // flag to run locally on AliAOD.root + AliAOD.VertexingHF.root
TString pathToLocalAODfiles="./AODfiles/LHC18r";         // path to find AOD files when running locally
bool runGridTest=false;                                  // flag to run a grid test: true (+runLocal=false). To run job on GRID: runGridTest=false, runLocal=false

// Alien output directory
TString gridWorkingDir="testNtupleCreatorDsPbPb2018";
TString gridOutputDir="output";

//Task configuration
TString cutFile="../cutobjects/DsDplusCuts_treecreator_PbPb2018_010_kINT7_kCentral.root"; 
TString cutFileDsTask="../cutobjects/DstoKKpiCuts_010_filttreecreator_Raa_kINT7_kCentral.root";
//______________________________________________

void RunAnalysisTreeCreator(TString configfilename, TString runMode)
{
    TGrid::Connect("alien://");
    // set if you want to run the analysis locally (true), or on grid (false)
    bool local = runLocal;
    // if you run on grid, specify test mode (true) or full grid model (false)
    bool gridTest = runGridTest;

    //load config
    YAML::Node config = YAML::LoadFile(configfilename.Data());
    string aliPhysVersion = config["AliPhysicVersion"].as<string>();
    string gridDataDir = config["datadir"].as<string>();
    string gridDataPattern = config["datapattern"].as<string>();
    string sSystem = config["system"].as<string>();
    int System=-1;
    if(sSystem=="pp")
        System = kpp;
    else if(sSystem=="pPb")
        System = kpPb;
    else if(sSystem=="PbPb")
        System = kPbPb;
    else {
        cerr << "Only pp, pPb and PbPb are supported, please check your config file. Exit" << endl;
    }
    vector<int> runlist = config["runs"].as<vector<int>>();
    const int nruns = runlist.size();
    bool isRunOnMC = strstr(gridDataDir.c_str(),"sim");

    // since we will compile a class, tell root where to look for headers
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    AliPhysicsSelectionTask *physseltask = reinterpret_cast<AliPhysicsSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s (%d)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"),isRunOnMC)));
    AliAnalysisTaskPIDResponse *pidResp = reinterpret_cast<AliAnalysisTaskPIDResponse *>(gInterpreter->ProcessLine(Form(".x %s (%d)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"),isRunOnMC)));

    if(System==kpPb || System==kPbPb) {
        AliMultSelectionTask *multSel = reinterpret_cast<AliMultSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"))));
        
        
        if(strstr(gridDataDir.c_str(),"LHC15o") || strstr(gridDataDir.c_str(),"LHC16i2")) multSel->SetAlternateOADBforEstimators("LHC15o-DefaultMC-HIJING");
        else if((strstr(gridDataDir.c_str(),"LHC18q") || strstr(gridDataDir.c_str(),"LHC18r") || strstr(gridDataDir.c_str(),"LHC19c")) && isRunOnMC) multSel->SetAlternateOADBFullManualBypassMC("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/data/OADB-LHC18q-DefaultMC-HIJING.root");
    }

    AliAnalysisTaskSEHFTreeCreator *task = reinterpret_cast<AliAnalysisTaskSEHFTreeCreator*>(gInterpreter->ProcessLine(Form(".x %s (%d,%d,\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d)",gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/treeHF/macros/AddTaskHFTreeCreator.C"),
    isRunOnMC, 1, "HFTreeCreator", cutFile.Data(),1,isRunOnMC,isRunOnMC,0,1,0,0,0,0,0,0,0,"AliHFTreeHandler::kNsigmaDetAndCombPID","AliHFTreeHandler::kNsigmaDetAndCombPID","AliHFTreeHandler::kNsigmaDetAndCombPID","AliHFTreeHandler::kNsigmaDetAndCombPID","AliHFTreeHandler::kNsigmaDetAndCombPID",
    "AliHFTreeHandler::kNsigmaDetAndCombPID","AliHFTreeHandler::kNsigmaDetAndCombPID","AliHFTreeHandler::kNsigmaDetAndCombPID","AliHFTreeHandler::kNsigmaDetAndCombPID","AliHFTreeHandler::kRedSingleTrackVars", false, false)));
    task->SetNsigmaTPCDataDrivenCorrection(0);

    AliAnalysisTaskSEDs *taskDs = reinterpret_cast<AliAnalysisTaskSEDs*>(gInterpreter->ProcessLine(Form(".x %s(%d,%d,%d,%d,%d,\"%s\",\"%s\")", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskDs.C"),1, 0, true,isRunOnMC,isRunOnMC,cutFileDsTask.Data(),"FilteringCuts")));
  
    if(System==kPbPb) {
        AliAnalysisTaskSECleanupVertexingHF *taskclean =reinterpret_cast<AliAnalysisTaskSECleanupVertexingHF *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskCleanupVertexingHF.C"))));
    }

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {

        // if you want to run locally, we need to define some input
        TChain* chainAOD = new TChain("aodTree");
        TChain *chainAODfriend = new TChain("aodTree");

        // add a few files to the chain (change this so that your local files are added)
        chainAOD->Add(Form("%s/AliAOD.root",pathToLocalAODfiles.Data()));
        chainAODfriend->Add(Form("%s/AliAOD.VertexingHF.root",pathToLocalAODfiles.Data()));

        chainAOD->AddFriend(chainAODfriend);


        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chainAOD);

    } 
    else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();

        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

        // make sure your source files get copied to grid
        // alienHandler->SetAdditionalLibs("custom.cxx custom.h");
        // alienHandler->SetAnalysisSource("custom.cxx");

        // select the aliphysics version. 
        alienHandler->SetAliPhysicsVersion(aliPhysVersion.data());

        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");

        // select the input data
        alienHandler->SetGridDataDir(gridDataDir.data());
        alienHandler->SetDataPattern(Form("%s/*AliAOD.root",gridDataPattern.data()));
        alienHandler->SetFriendChainName("AliAOD.VertexingHF.root");

        // MC has no prefix, data has prefix 000
        if(!isRunOnMC)
            alienHandler->SetRunPrefix("000");

        for(int k=0; k<nruns; k++)
            alienHandler->AddRunNumber(runlist[k]);
        alienHandler->SetNrunsPerMaster(nruns);

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(1);
        alienHandler->SetExecutable("myTask.sh");

        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(true);
        alienHandler->SetKeepLogs(true);

        // merging: run with true to merge on grid
        // after re-running the jobs in SetRunMode("terminate")
        // (see below) mode, set SetMergeViaJDL(false)
        // to collect final results
        alienHandler->SetMaxMergeStages(3); //2, 3
        alienHandler->SetMergeViaJDL(true);

        // define the output folders
        alienHandler->SetGridWorkingDir(gridWorkingDir.Data());
        alienHandler->SetGridOutputDir(gridOutputDir.Data());

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);

        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        }
        else {
            // else launch the full grid analysis
            alienHandler->SetRunMode(runMode.Data()); //terminate
            mgr->StartAnalysis("grid");
        }
    }
}
