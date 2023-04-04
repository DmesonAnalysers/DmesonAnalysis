'''
Script for validation of tree creator output
run: python ValidateTreeCreator inputfile.root inputdir inputlist
if the output is from MC, use --mc option
'''

import sys
import argparse
import pandas as pd
import numba
import uproot
from ROOT import TFile, TCanvas, TH1F, TLegend # pylint: disable=import-error, no-name-in-module
from ROOT import kBlack, kRed, kFullCircle, kOpenSquare # pylint: disable=import-error, no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error


@numba.njit
def FidAccSel(arrayPt, arrayY):
    '''
    method for the appliction of the fiducial acceptance selection
    '''
    arrayIsSel = []
    for icand, pt in enumerate(arrayPt):
        if pt > 5:
            if abs(arrayY[icand]) < 0.8:
                arrayIsSel.append(True)
            else:
                arrayIsSel.append(False)
        else:
            yfid = -0.2/15 * pt**2 + 1.9/15 * pt + 0.5
            if abs(arrayY[icand]) < yfid:
                arrayIsSel.append(True)
            else:
                arrayIsSel.append(False)
    return arrayIsSel


# Main function
Vars = {'inv_mass':0}

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('--mc', action='store_true', help='flag for MC outputs')
parser.add_argument('inputfile', metavar='text', help='input root file name')
parser.add_argument('inputdir', metavar='text', help='task directory name')
parser.add_argument('inputlist', metavar='text', help='task list name')
args = parser.parse_args()

if args.mc:
    print('MC validation not yet implemented. Exit')
    sys.exit()

SetGlobalStyle(padleftmargin=0.15, padtopmargin=0.05, titlesize=0.045, labelsize=0.04)

hMassTask, hMassTreeCreator, cComparison = ({} for iDic in range(3))

infile = TFile.Open(args.inputfile)
listTask = infile.Get('{0}/{1}'.format(args.inputdir, args.inputlist))
hEv = listTask.FindObject('hNEvents')
sparse = listTask.FindObject('fnSparse')
for iVar in Vars:
    hMassTask[iVar] = sparse.Projection(Vars[iVar])
    SetObjectStyle(hMassTask[iVar], color=kBlack, markerstyle=kFullCircle)

dfD = uproot.open(args.inputfile)['PWGHF_TreeCreator/tree_Ds'].arrays(library='pd') #tree_event_char
dfEv = uproot.open(args.inputfile)['PWGHF_TreeCreator/tree_event_char'].arrays(library='pd')
dfEvSel = dfEv.query('is_ev_rej==0')

dfMerged = pd.merge(dfD, dfEvSel, on=['run_number', 'ev_id'])

isSel = FidAccSel(dfMerged['pt_cand'].to_numpy(), dfMerged.y_cand.to_numpy())
dfMerged['fidacc'] = isSel
dfMergedSel = dfMerged.query('fidacc==True')

leg = TLegend(0.4, 0.6, 0.8, 0.8)
for counter, iVar in enumerate(Vars):
    nbins = hMassTask[iVar].GetNbinsX()
    mmin = hMassTask[iVar].GetBinLowEdge(1)
    mmax = hMassTask[iVar].GetBinLowEdge(nbins)+hMassTask[iVar].GetBinWidth(nbins)
    hMassTreeCreator[iVar] = TH1F('hMassTreeCreator{0}'.format(iVar), '', nbins, mmin, mmax)
    for mass in dfMergedSel[iVar].to_numpy():
        hMassTreeCreator[iVar].Fill(mass)
    SetObjectStyle(hMassTreeCreator[iVar], color=kRed, markerstyle=kOpenSquare)
    cComparison[iVar] = TCanvas('cComparison{0}'.format(iVar), '', 800, 800)
    ymax = hMassTask[iVar].GetMaximum()*2
    cComparison[iVar].DrawFrame(mmin, 0.1, mmax, ymax, ';{0}; Counts'.format(iVar))
    hMassTask[iVar].Draw('Esame')
    hMassTreeCreator[iVar].Draw('Esame')
    if counter == 0:
        leg.AddEntry(hMassTask[iVar], 'task', 'p')
        leg.AddEntry(hMassTreeCreator[iVar], 'tree creator', 'p')
    leg.Draw()
    cComparison[iVar].SaveAs('{0}.pdf'.format(iVar))

hSelEv = TH1F('hEvSel', ';;selected events', 2, 0.5, 2.5)
hSelEv.GetXaxis().SetBinLabel(1, 'task')
hSelEv.GetXaxis().SetBinLabel(2, 'tree creator')
hSelEv.SetBinContent(1, hEv.GetBinContent(5))
hSelEv.SetBinContent(2, len(dfEvSel))
SetObjectStyle(hSelEv, color=kBlack, markerstyle=kFullCircle)

cEvents = TCanvas('cEvents', '', 800, 800)
hSelEv.Draw('E')
cEvents.SaveAs('selected_events.pdf')

input('Press enter to exit')
