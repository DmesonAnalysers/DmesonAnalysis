AliAnalysisTaskCheckDprotonPairGen *AddTaskCheckDprotonPairGen()
{
    // Creates, configures and attaches to the train the task for QA of ITS standalone tracks
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AliAnalysisTaskCheckDprotonPairGen", "No analysis manager to connect to.");
        return NULL;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler())
    {
        ::Error("AliAnalysisTaskCheckDprotonPairGen", "This task requires an input event handler");
        return NULL;
    }

    // Create and configure the task
    AliAnalysisTaskCheckDprotonPairGen *task = new AliAnalysisTaskCheckDprotonPairGen();
    mgr->AddTask(task);

    // Create ONLY the output containers for the data produced by the task.
    // Get and connect other common input/output containers via the manager as below
    //==============================================================================
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    outputFileName += ":DprotonPairGen";

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutputDprotonPairGen",
                                                              TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              outputFileName);

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput1);
    return task;
}
