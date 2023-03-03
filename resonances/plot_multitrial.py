'''
python script to plot the results from multitrial systematics
run: python plot_multitrial.py
'''
import sys
import argparse
import yaml
import numpy as np
from ROOT import TFile, TH2F, TH1F, TCanvas, TLatex, kBlack, kRed, kAzure, kGray, TLine
sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
SetGlobalStyle(padtopmargin=0.05, padleftmargin=0.18,padbottommargin=0.15,
               palette=55, labelsize=0.04, titlesize=0.06,
               labeloffset=0.008, titleoffsety=1.55, titleoffsetx=1.2, titleoffsetz=1.,
               opttitle=0, optstat=0)

reso = 435
mult = 'MB'
outlabel = '_limitedges_wchebpol2'

#______________________________________________________________________________
# Load data
if reso == 435:
    infile_syst = f'/home/stefano/Desktop/Analyses/DmesonAnalysis/resonances/Multitrial_Ds2starplus_pt2.0-24.0_{mult}{outlabel}.root'
    infile_ref = f'/home/stefano/Desktop/cernbox/Ds_resonances/raw_yield/mass_Ds2starplus_pt2.0-24.0_{mult}.root'
    resoname = 'Ds2starplus'
    reso_label = 'D_{s2}^{*+} #rightarrow D^{+} K_{S}^{0}'
else:
    infile_syst = f'/home/stefano/Desktop/Analyses/DmesonAnalysis/resonances/Multitrial_Ds1plus_pt2.0-24.0_{mult}{outlabel}.root'
    infile_ref = f'/home/stefano/Desktop/cernbox/Ds_resonances/raw_yield/mass_Ds1plus_pt2.0-24.0_{mult}.root'
    resoname = 'Ds1plus'
    reso_label = 'D_{s1}^{+} #rightarrow D^{*+} K_{S}^{0}'
print(f'\033[1m\033[92m Loading multitrial file {infile_syst} \033[0m')
print(f'\033[1m\033[92m Loading reference file {infile_ref} \033[0m')

#_____________________________________________________________________________
# Multitrial output
infile_root = TFile.Open(infile_syst)
hraw = infile_root.Get('h_rawyields')
hgamma = infile_root.Get('h_gamma')
hsignif = infile_root.Get('h_signif')
hsigma = infile_root.Get('h_sigmas')
hchi2 = infile_root.Get('h_chi2')
hsyst = infile_root.Get('h_syst')
ntrials = hsignif.GetXaxis().GetXmax()

for h in [hraw, hgamma, hsignif, hsigma, hchi2, hsyst]:
    col = kAzure+4 if ('syst' in h.GetName() or 'raw' in h.GetName()) else kBlack
    fillstyle = 3345 if 'syst' in h.GetName() else 0
    SetObjectStyle(h, color=col, markerstyle=20, markersize=1.,
                  linecolor=col, markercolor=col, fillalpha=0.2, fillstyle=fillstyle,
                  linewidth=1)

#_____________________________________________________________________________
# Reference raw yield and thresholds
ref_file = TFile.Open(infile_ref)
hraw_ref = ref_file.Get('h_rawyields')
hsigma_ref = ref_file.Get('h_sigmas')

rawyield_ref = hraw_ref.GetBinContent(1)
sigma_ref = hsigma_ref.GetBinContent(1)
line_rawref = TLine(0, rawyield_ref, ntrials, rawyield_ref)
line_sigmaref = TLine(0, sigma_ref, ntrials, sigma_ref)
line_syst = TLine(rawyield_ref, 0, rawyield_ref, hsyst.GetMaximum()*1.05)
line_chi2 = TLine(0, 2., ntrials, 2.) #  chi2/ndf = 2 threshold
l_trh = TLine(0., 3., ntrials, 3.) #  signif = 3 threshold

line_rawref.SetLineColor(kRed+1)
line_rawref.SetLineStyle(1)
line_rawref.SetLineWidth(3)
line_sigmaref.SetLineColor(kRed+1)
line_sigmaref.SetLineStyle(1)
line_sigmaref.SetLineWidth(3)
line_syst.SetLineColor(kRed+1)
line_syst.SetLineStyle(1)
line_syst.SetLineWidth(3)
line_chi2.SetLineColor(kGray+2)
line_chi2.SetLineStyle(2)
l_trh.SetLineColor(kGray+2)
l_trh.SetLineStyle(2)
l_trh.SetLineWidth(2)

canvas = TCanvas('canvas', 'canvas', 1000, 1800)
canvas.Divide(2, 3)

hist_to_plot = [hraw, hsyst, hgamma, hsignif, hsigma, hchi2]
lines = [line_rawref, line_syst, None, l_trh, line_sigmaref, line_chi2]
titles = ['; trial; raw yield',
          '; raw yield; entries',
          '; trial; #Gamma (GeV/#it{c}^{2})',
          '; trial; significance',
          '; trial; #sigma (GeV/#it{c}^{2})',
          '; trial; #chi^{2}/ndf']
for pad in range(1, 7):
    hist = hist_to_plot[pad-1]

    if pad == 2:
        xcenter = hist.GetXaxis().GetBinCenter(hist.GetMaximumBin())
        xmin = xcenter - 20*hist.GetRMS()
        xmax = xcenter + 20*hist.GetRMS()
    else:
        xmin = 0.
        xmax = ntrials
    if pad in [1, 2, 3, 5]:
        ymin = hist.GetMinimum()*0.5
        ymax = hist.GetMaximum()*1.8
    else:
        ymin = 0.
        ymax = 5 if pad == 6 else 10
    hFrame = canvas.cd(pad).DrawFrame(xmin, ymin, xmax, ymax, titles[pad-1])
    hFrame.GetXaxis().SetNdivisions(505)

    draw_opt = 'pz same' if pad != 2 else 'hist same'
    hist.Draw(draw_opt)
    if lines[pad-1]:
        lines[pad-1].Draw('same')

    if pad == 2:
        rms = hsyst.GetRMS()
        mean = hsyst.GetMean()
        max_bin = 10000
        for i in range(1, hsyst.GetNbinsX()+1):
            if hsyst.GetBinContent(i) > 0.:
                max_bin = i
        min_bin = 0
        for i in range(1, hsyst.GetNbinsX()+1):
            if hsyst.GetBinContent(i) > 0.:
                min_bin = i
                break

        shift = (hsyst.GetBinCenter(max_bin) - hsyst.GetBinCenter(min_bin))/np.sqrt(12)
        syst = np.sqrt((rms/mean)**2 + (shift/mean)**2)
        lat_syst = TLatex()
        lat_syst.SetTextSize(0.04)
        lat_syst.SetNDC()
        lat_syst.SetTextFont(42)
        lat_syst.SetTextColor(kBlack)
        lat_syst.DrawLatexNDC(0.55, 0.88, f'mean = {mean:.2f}')
        lat_syst.DrawLatexNDC(0.55, 0.81, f'RMS = {rms:.2f} ({rms*100/mean:.0f}%)')
        lat_syst.DrawLatexNDC(0.55, 0.73,
                              f'#frac{{max-min}}{{#sqrt{{12}}}} = {shift:.2f} ({shift*100/mean:.0f}%)')
        lat_syst.DrawLatexNDC(0.2, 0.88, f'{reso_label} ({mult})')

canvas.Update()
outfile = TFile.Open(f'./multitrial_syst_{resoname}_{mult}{outlabel}.root', 'recreate')
for h in hist_to_plot:
    h.Write()
canvas.Write()
canvas.SaveAs(f'./multitrial_syst_{resoname}_{mult}{outlabel}.pdf')
outfile.Close()
input('Press enter to exit')
