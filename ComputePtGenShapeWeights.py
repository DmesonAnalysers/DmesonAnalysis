'''
Script for the computation of pT shape weights
'''

import array
from ROOT import TFile, TSpline3  # pylint: disable=import-error,no-name-in-module
from ReadModel import ReadTAMU, ReadPHSD, ReadGossiaux, ReadCatania

infileGenPtShape = TFile.Open('ptweights/GenPtShape_LHC19c3a.root')
hPtGen = infileGenPtShape.Get('hPtPromptStep0')
hPtGen.SetName('hPtPythia')
hPtGen.SetDirectory(0)
nptbins = hPtGen.GetNbinsX()
ptlims = array.array('d')
for iPt in range(nptbins):
    ptlims.append(hPtGen.GetBinLowEdge(iPt+1))
ptlims.append(hPtGen.GetBinLowEdge(nptbins)+hPtGen.GetBinWidth(nptbins))

infileFONLL = TFile.Open(
    'models/D0DplusDstarPredictions_502TeV_y05_all_021016_BDShapeCorrected.root')
hPtFONLLcent = infileFONLL.Get('hDsPhipitoKkpipred_central')
hPtFONLLmin = infileFONLL.Get('hDsPhipitoKkpipred_min')
hPtFONLLmax = infileFONLL.Get('hDsPhipitoKkpipred_max')

TAMU = ReadTAMU('models/Ds_TAMU_RAA_5TeV_010.txt')
PHSD = ReadPHSD('models/Ds_PHSD_RAA_5TeV_010.txt')
Gossiaux = ReadGossiaux('models/D_Gossiaux_RAA_5TeV_010.txt')
Catania = ReadCatania('models/Ds_Catania_RAA_5TeV_010.dat')

sTAMU, sGossiaux, sPHSD, sCatania = ({} for iDic in range(4))

for iRaaTAMU in TAMU:
    if iRaaTAMU != 'PtCent':
        sTAMU[iRaaTAMU] = TSpline3('sTAMU%s' % iRaaTAMU, array.array(
            'd', TAMU['PtCent']), array.array('d', TAMU[iRaaTAMU]), len(TAMU['PtCent']))
for iRaaPHSD in PHSD:
    if iRaaPHSD != 'PtCent':
        sPHSD[iRaaPHSD] = TSpline3('sPHSD%s' % iRaaPHSD, array.array(
            'd', PHSD['PtCent']), array.array('d', PHSD[iRaaPHSD]), len(PHSD['PtCent']))
for iRaaGossiaux in Gossiaux:
    if iRaaGossiaux != 'PtCent':
        sGossiaux[iRaaGossiaux] = TSpline3('sGossiaux%s' % iRaaGossiaux, array.array(
            'd', Gossiaux['PtCent']), array.array('d', Gossiaux[iRaaGossiaux]), len(Gossiaux['PtCent']))
for iRaaCatania in Catania:
    if iRaaCatania != 'PtCent':
        sCatania[iRaaCatania] = TSpline3('sCatania%s' % iRaaCatania, array.array(
            'd', Catania['PtCent']), array.array('d', Catania[iRaaCatania]), len(Catania['PtCent']))

hPtFONLLcent = hPtFONLLcent.Rebin(nptbins, 'hPtFONLLcent', ptlims)
hPtFONLLmin = hPtFONLLmin.Rebin(nptbins, 'hPtFONLLmin', ptlims)
hPtFONLLmax = hPtFONLLmax.Rebin(nptbins, 'hPtFONLLmax', ptlims)

hPtFONLLtimesTAMUcent = hPtFONLLcent.Clone('hPtFONLLtimesTAMUcent')
hPtFONLLtimesTAMUmin = hPtFONLLmin.Clone('hPtFONLLtimesTAMUmin')
hPtFONLLtimesTAMUmax = hPtFONLLmax.Clone('hPtFONLLtimesTAMUmax')

hPtFONLLtimesPHSDcent = hPtFONLLcent.Clone('hPtFONLLtimesPHSDcent')
hPtFONLLtimesPHSDmin = hPtFONLLmin.Clone('hPtFONLLtimesPHSDmin')
hPtFONLLtimesPHSDmax = hPtFONLLmax.Clone('hPtFONLLtimesPHSDmax')

hPtFONLLtimesGossiauxcent = hPtFONLLcent.Clone('hPtFONLLtimesGossiauxcent')
hPtFONLLtimesGossiauxmin = hPtFONLLmin.Clone('hPtFONLLtimesGossiauxmin')
hPtFONLLtimesGossiauxmax = hPtFONLLmax.Clone('hPtFONLLtimesGossiauxmax')

hPtFONLLtimesCataniacent = hPtFONLLcent.Clone('hPtFONLLtimesCataniacent')
hPtFONLLtimesCataniamin = hPtFONLLmin.Clone('hPtFONLLtimesCataniamin')
hPtFONLLtimesCataniamax = hPtFONLLmax.Clone('hPtFONLLtimesCataniamax')

for iPt in range(nptbins):
    pt = hPtFONLLcent.GetBinCenter(iPt+1)
    hPtFONLLtimesTAMUcent.SetBinContent(
        iPt+1, hPtFONLLcent.GetBinContent(iPt+1)*(sTAMU['Max'].Eval(pt)+sTAMU['Min'].Eval(pt))/2)
    hPtFONLLtimesTAMUmin.SetBinContent(
        iPt+1, hPtFONLLmin.GetBinContent(iPt+1)*(sTAMU['Max'].Eval(pt)+sTAMU['Min'].Eval(pt))/2)
    hPtFONLLtimesTAMUmax.SetBinContent(
        iPt+1, hPtFONLLmax.GetBinContent(iPt+1)*(sTAMU['Max'].Eval(pt)+sTAMU['Min'].Eval(pt))/2)

    hPtFONLLtimesPHSDcent.SetBinContent(
        iPt+1, hPtFONLLcent.GetBinContent(iPt+1)*sPHSD['Cent'].Eval(pt))
    hPtFONLLtimesPHSDmin.SetBinContent(
        iPt+1, hPtFONLLmin.GetBinContent(iPt+1)*sPHSD['Cent'].Eval(pt))
    hPtFONLLtimesPHSDmax.SetBinContent(
        iPt+1, hPtFONLLmax.GetBinContent(iPt+1)*sPHSD['Cent'].Eval(pt))

    hPtFONLLtimesGossiauxcent.SetBinContent(
        iPt+1, hPtFONLLcent.GetBinContent(iPt+1)*sGossiaux['ColRad'].Eval(pt))
    hPtFONLLtimesGossiauxmin.SetBinContent(
        iPt+1, hPtFONLLmin.GetBinContent(iPt+1)*sGossiaux['ColRad'].Eval(pt))
    hPtFONLLtimesGossiauxmax.SetBinContent(
        iPt+1, hPtFONLLmax.GetBinContent(iPt+1)*sGossiaux['ColRad'].Eval(pt))

    hPtFONLLtimesCataniacent.SetBinContent(
        iPt+1, hPtFONLLcent.GetBinContent(iPt+1)*sCatania['Cent'].Eval(pt))
    hPtFONLLtimesCataniamin.SetBinContent(
        iPt+1, hPtFONLLmin.GetBinContent(iPt+1)*sCatania['Cent'].Eval(pt))
    hPtFONLLtimesCataniamax.SetBinContent(
        iPt+1, hPtFONLLmax.GetBinContent(iPt+1)*sCatania['Cent'].Eval(pt))

hPtGen.Scale(1./hPtGen.Integral())

hPtFONLLcent.Scale(1./hPtFONLLcent.Integral())
hPtFONLLmin.Scale(1./hPtFONLLmin.Integral())
hPtFONLLmax.Scale(1./hPtFONLLmax.Integral())

hPtFONLLtimesTAMUcent.Scale(1./hPtFONLLtimesTAMUcent.Integral())
hPtFONLLtimesTAMUmin.Scale(1./hPtFONLLtimesTAMUmin.Integral())
hPtFONLLtimesTAMUmax.Scale(1./hPtFONLLtimesTAMUmax.Integral())

hPtFONLLtimesPHSDcent.Scale(1./hPtFONLLtimesPHSDcent.Integral())
hPtFONLLtimesPHSDmin.Scale(1./hPtFONLLtimesPHSDmin.Integral())
hPtFONLLtimesPHSDmax.Scale(1./hPtFONLLtimesPHSDmax.Integral())

hPtFONLLtimesGossiauxcent.Scale(1./hPtFONLLtimesGossiauxcent.Integral())
hPtFONLLtimesGossiauxmin.Scale(1./hPtFONLLtimesGossiauxmin.Integral())
hPtFONLLtimesGossiauxmax.Scale(1./hPtFONLLtimesGossiauxmax.Integral())

hPtFONLLtimesCataniacent.Scale(1./hPtFONLLtimesCataniacent.Integral())
hPtFONLLtimesCataniamin.Scale(1./hPtFONLLtimesCataniamin.Integral())
hPtFONLLtimesCataniamax.Scale(1./hPtFONLLtimesCataniamax.Integral())

hPtWeightsFONLLcent = hPtFONLLcent.Clone('hPtWeightsFONLLcent')
hPtWeightsFONLLmin = hPtFONLLmin.Clone('hPtWeightsFONLLmin')
hPtWeightsFONLLmax = hPtFONLLmax.Clone('hPtWeightsFONLLmax')
hPtWeightsFONLLcent.Divide(hPtFONLLcent, hPtGen)
hPtWeightsFONLLmin.Divide(hPtFONLLmin, hPtGen)
hPtWeightsFONLLmax.Divide(hPtFONLLmax, hPtGen)

hPtWeightsFONLLtimesTAMUcent = hPtFONLLcent.Clone('hPtWeightsFONLLtimesTAMUcent')
hPtWeightsFONLLtimesTAMUmin = hPtFONLLmin.Clone('hPtWeightsFONLLtimesTAMUmin')
hPtWeightsFONLLtimesTAMUmax = hPtFONLLmax.Clone('hPtWeightsFONLLtimesTAMUmax')
hPtWeightsFONLLtimesTAMUcent.Divide(hPtFONLLtimesTAMUcent, hPtGen)
hPtWeightsFONLLtimesTAMUmin.Divide(hPtFONLLtimesTAMUmin, hPtGen)
hPtWeightsFONLLtimesTAMUmax.Divide(hPtFONLLtimesTAMUmax, hPtGen)

hPtWeightsFONLLtimesPHSDcent = hPtFONLLcent.Clone('hPtWeightsFONLLtimesPHSDcent')
hPtWeightsFONLLtimesPHSDmin = hPtFONLLmin.Clone('hPtWeightsFONLLtimesPHSDmin')
hPtWeightsFONLLtimesPHSDmax = hPtFONLLmax.Clone('hPtWeightsFONLLtimesPHSDmax')
hPtWeightsFONLLtimesPHSDcent.Divide(hPtFONLLtimesPHSDcent, hPtGen)
hPtWeightsFONLLtimesPHSDmin.Divide(hPtFONLLtimesPHSDmin, hPtGen)
hPtWeightsFONLLtimesPHSDmax.Divide(hPtFONLLtimesPHSDmax, hPtGen)

hPtWeightsFONLLtimesGossiauxcent = hPtFONLLcent.Clone('hPtWeightsFONLLtimesGossiauxcent')
hPtWeightsFONLLtimesGossiauxmin = hPtFONLLmin.Clone('hPtWeightsFONLLtimesGossiauxmin')
hPtWeightsFONLLtimesGossiauxmax = hPtFONLLmax.Clone('hPtWeightsFONLLtimesGossiauxmax')
hPtWeightsFONLLtimesGossiauxcent.Divide(hPtFONLLtimesGossiauxcent, hPtGen)
hPtWeightsFONLLtimesGossiauxmin.Divide(hPtFONLLtimesGossiauxmin, hPtGen)
hPtWeightsFONLLtimesGossiauxmax.Divide(hPtFONLLtimesGossiauxmax, hPtGen)

hPtWeightsFONLLtimesCataniacent = hPtFONLLcent.Clone('hPtWeightsFONLLtimesCataniacent')
hPtWeightsFONLLtimesCataniamin = hPtFONLLmin.Clone('hPtWeightsFONLLtimesCataniamin')
hPtWeightsFONLLtimesCataniamax = hPtFONLLmax.Clone('hPtWeightsFONLLtimesCataniamax')
hPtWeightsFONLLtimesCataniacent.Divide(hPtFONLLtimesCataniacent, hPtGen)
hPtWeightsFONLLtimesCataniamin.Divide(hPtFONLLtimesCataniamin, hPtGen)
hPtWeightsFONLLtimesCataniamax.Divide(hPtFONLLtimesCataniamax, hPtGen)

outfile = TFile('ptweights/PtWeigths_LHC19c3a.root', 'recreate')
# spectrum shapes
hPtGen.Write()
hPtFONLLcent.Write()
hPtFONLLmin.Write()
hPtFONLLmax.Write()
hPtFONLLtimesTAMUcent.Write()
hPtFONLLtimesTAMUmin.Write()
hPtFONLLtimesTAMUmax.Write()
hPtFONLLtimesPHSDcent.Write()
hPtFONLLtimesPHSDmin.Write()
hPtFONLLtimesPHSDmax.Write()
hPtFONLLtimesGossiauxcent.Write()
hPtFONLLtimesGossiauxmin.Write()
hPtFONLLtimesGossiauxmax.Write()
hPtFONLLtimesCataniacent.Write()
hPtFONLLtimesCataniamin.Write()
hPtFONLLtimesCataniamax.Write()
# weights
hPtWeightsFONLLcent.Write()
hPtWeightsFONLLmin.Write()
hPtWeightsFONLLmax.Write()
hPtWeightsFONLLtimesTAMUcent.Write()
hPtWeightsFONLLtimesTAMUmin.Write()
hPtWeightsFONLLtimesTAMUmax.Write()
hPtWeightsFONLLtimesPHSDcent.Write()
hPtWeightsFONLLtimesPHSDmin.Write()
hPtWeightsFONLLtimesPHSDmax.Write()
hPtWeightsFONLLtimesGossiauxcent.Write()
hPtWeightsFONLLtimesGossiauxmin.Write()
hPtWeightsFONLLtimesGossiauxmax.Write()
hPtWeightsFONLLtimesCataniacent.Write()
hPtWeightsFONLLtimesCataniamin.Write()
hPtWeightsFONLLtimesCataniamax.Write()
outfile.Close()
