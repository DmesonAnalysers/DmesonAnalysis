# DmesonCutSelectionPbPb

#run significance optimisation:
#python ScanSignificanceSparse.py configfiles/cfgFileName.yml outFileName.root

#project significance ntuple (setting a minimum of the accepted significance and efficiency):
#python ProjectSignifNtuple.py configfiles/cfgFileName.yml inFileName.root PtMin PtMax minSignificance minEffPrompt

#produce plots with chosen cuts:
#python ProjectDplusDsSparse.py configfiles/cfgFileName.yml configfiles/cutSetFile.yml outFileName