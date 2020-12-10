'''
Script for HEP data reading
'''

import math
from ROOT import TFile, TGraphAsymmErrors # pylint: disable=import-error,no-name-in-module

def ReadHepDataROOT(HepDataFileName, tablenum):
    infile = TFile(HepDataFileName)
    table = infile.Get('Table %d' % tablenum)
    Hist1D_y1 = table.Get('Hist1D_y1')
    Hist1D_y1_e1 = table.Get('Hist1D_y1_e1')
    Hist1D_y1_e2plus = table.Get('Hist1D_y1_e2plus')
    Hist1D_y1_e2minus = table.Get('Hist1D_y1_e2minus')
    Hist1D_y1_e3 = table.Get('Hist1D_y1_e3')
    Hist1D_y1.SetDirectory(0)
    Hist1D_y1_e1.SetDirectory(0)
    Hist1D_y1_e2plus.SetDirectory(0)
    Hist1D_y1_e2minus.SetDirectory(0)
    Hist1D_y1_e3.SetDirectory(0)
    infile.Close()

    graphSyst = TGraphAsymmErrors()
    histoStat = Hist1D_y1.Clone('histoStat')

    for iPt in range(Hist1D_y1.GetNbinsX()):
        histoStat.SetBinContent(iPt+1, Hist1D_y1.GetBinContent(iPt+1))
        histoStat.SetBinError(iPt+1, Hist1D_y1_e1.GetBinContent(iPt+1))
        systlow = math.sqrt(Hist1D_y1_e2minus.GetBinContent(iPt+1)**2+Hist1D_y1_e3.GetBinContent(iPt+1)**2)
        systhigh = math.sqrt(Hist1D_y1_e2plus.GetBinContent(iPt+1)**2+Hist1D_y1_e3.GetBinContent(iPt+1)**2)
        graphSyst.SetPoint(iPt, histoStat.GetBinCenter(iPt+1), histoStat.GetBinContent(iPt+1))
        graphSyst.SetPointError(iPt, histoStat.GetBinWidth(iPt+1)/2, histoStat.GetBinWidth(iPt+1)/2, systlow, systhigh)

    return histoStat, graphSyst
