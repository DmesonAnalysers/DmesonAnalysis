#*******************************************************************************************#
# python script for the projection of D+ and Ds+ mesons THnSparses                          #
# run: python ProjectDplusDsSparse.py cfgFileName.yml cutSetFileName.yml outFileName.root   #
# author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino                            #
#*******************************************************************************************#

from ROOT import TFile, TCanvas, TH1F
from ROOT import gROOT
import yaml
import sys

cfgFileName = sys.argv[1]
cutSetFileName = sys.argv[2]
outFileName = sys.argv[3]

with open(cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile)

isMC = inputCfg['isMC']
infileData = TFile(inputCfg['filename'])
indirData = infileData.Get(inputCfg['dirname'])
inlistData = indirData.Get(inputCfg['listname'])
sMassPtCutVars = inlistData.FindObject(inputCfg['sparsenameAll'])
if isMC:
    sMassPtCutVarsPrompt = inlistData.FindObject(inputCfg['sparsenamePrompt'])
    sMassPtCutVarsFD = inlistData.FindObject(inputCfg['sparsenameFD'])
    sGenPrompt = inlistData.FindObject(inputCfg['sparsenameGenPrompt'])
    sGenFD = inlistData.FindObject(inputCfg['sparsenameGenFD'])
normCounter = indirData.Get(inputCfg['normname'])
hEv = inlistData.FindObject(inputCfg['histoevname'])

with open(cutSetFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile)

cutVars = cutSetCfg['cutvars']

outfile = TFile(outFileName,'RECREATE')

for iPt in range(0,len(cutVars['Pt']['min'])):
    for iVar in cutVars:
        if iVar == 'InvMass':
            continue
        binMin = sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['min'][iPt]*1.0001)
        binMax = sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['max'][iPt]*0.9999)
        sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin,binMax)
        if isMC:
            sMassPtCutVarsPrompt.GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin,binMax)
            sMassPtCutVarsFD.GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin,binMax)
    for iVar in cutVars:
        hVar = sMassPtCutVars.Projection(cutVars[iVar]['axisnum'])
        hVar.SetName('h%s_%0.f_%0.f' % ( cutVars[iVar]['name'], cutVars['Pt']['min'][iPt] , cutVars['Pt']['max'][iPt]))    
        outfile.cd()
        hVar.Write()
        if isMC:
            hVarPrompt = sMassPtCutVarsPrompt.Projection(cutVars[iVar]['axisnum'])
            hVarPrompt.SetName('hPrompt%s_%0.f_%0.f' % ( cutVars[iVar]['name'], cutVars['Pt']['min'][iPt] , cutVars['Pt']['max'][iPt]))    
            hVarPrompt.Write()
            hVarFD = sMassPtCutVarsFD.Projection(cutVars[iVar]['axisnum'])
            hVarFD.SetName('hFD%s_%0.f_%0.f' % ( cutVars[iVar]['name'], cutVars['Pt']['min'][iPt] , cutVars['Pt']['max'][iPt]))    
            hVarFD.Write()
    if isMC:
        binGenMin = sGenPrompt.GetAxis(0).FindBin(cutVars['Pt']['min'][iPt]*1.0001)
        binGenMax = sGenPrompt.GetAxis(0).FindBin(cutVars['Pt']['max'][iPt]*0.9999)
        sGenPrompt.GetAxis(0).SetRange(binGenMin,binGenMax)
        sGenFD.GetAxis(0).SetRange(binGenMin,binGenMax)
        hGenPtPrompt = sGenPrompt.Projection(0)
        hGenPtPrompt.SetName('hPrompGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt] , cutVars['Pt']['max'][iPt]))    
        hGenPtPrompt.Write()
        hGenPtFD = sGenFD.Projection(0)
        hGenPtFD.SetName('hFDGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt] , cutVars['Pt']['max'][iPt]))    
        hGenPtFD.Write()

hEvForNorm = TH1F("hEvForNorm",";;Number of events",2,0.,2.)
hEvForNorm.GetXaxis().SetBinLabel(1,"norm counter")
hEvForNorm.GetXaxis().SetBinLabel(2,"accepted events")
hEvForNorm.SetBinContent(1,normCounter.GetNEventsForNorm())
hEvForNorm.SetBinContent(2,hEv.GetBinContent(10))

hEvForNorm.Write()
outfile.Close()
