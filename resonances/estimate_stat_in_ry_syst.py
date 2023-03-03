'''
python script to generate the multitrial ditribution for systematic checks
run: python generate_multitrial.py
'''
import sys
import numpy as np
from ROOT import TFile, TH1F, TGraphErrors, kRed, TF1, kBlack, TCanvas, TLatex
sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
SetGlobalStyle(padtopmargin=0.05, padleftmargin=0.18,padbottommargin=0.15,
               palette=55, labelsize=0.04, titlesize=0.06,
               labeloffset=0.008, titleoffsety=1.55, titleoffsetx=1.2, titleoffsetz=1.,
               opttitle=0, optstat=0)
name_reso = "Ds2starplus"
trigger = "HM"
outlabel = ""

#______________________________________________________________________________
# Load reference file
entries = [5, 30, 60, 70, 80] # HM: 5, 30, 60, 70, 80; MB: 1, 5, 10, 20, 30
tfiles = []
gRMS = TGraphErrors()
SetObjectStyle(gRMS, color=kBlack, markerstyle=20, markersize=0.5)
htrials = []
rms, rms_unc = [], []
for ientry, entry in enumerate(entries):
    for trial in range(max(entries)):
        if trigger == "HM" and trial in [12, 62, 78, 79, 80]: # skip bad trials HM: [5, 30, 60, 70, 80]
            continue
        if trial <= entry:
            print(f"trial {trial} in range 0-{entry}")
            tfile = TFile(f'./Multitrial_{name_reso}_pt2.0-24.0_{trigger}_test_bash_htrial_{trial}.root', "READ")
            htrials.append(tfile.Get("h_syst"))
            htrials[-1].SetDirectory(0)
        else:
            continue

    if not htrials:
        continue
    else:
        print(f"Found {len(htrials)} histograms")
    hfinal = TH1F(f"hfinal_{entry}", "hfinal", 1000, 0, 1)
    for ibin in range(1, hfinal.GetNbinsX()+1):
        hfinal.SetBinContent(ibin, 0)
    hfinal.SetDirectory(0)
    SetObjectStyle(hfinal, color=kRed, markerstyle=20, markersize=0.5)
    for h in htrials:
        hfinal.Fill(h.GetRMS()/h.GetMean())
    rms.append(hfinal.GetMean())
    rms_unc.append(hfinal.GetMeanError()/np.sqrt(len(htrials)))
    htrials = []

f1 = TF1("f1", "[0] + [1]/sqrt(x)", 0, 100)
canvas = TCanvas("canvas", "canvas", 800, 800)
ymax = max(rms) + 0.05
ymin = min(rms) - 0.02
hFrame = canvas.DrawFrame(0, ymin, max(entries)+5, ymax, ";N_{trials};Average RMS/Mean")
for i, (entry, rms, rms_unc) in enumerate(zip(entries, rms, rms_unc)):
    gRMS.SetPoint(i, entry, rms)
    gRMS.SetPointError(i, 0, rms_unc)
gRMS.Draw('p same')
gRMS.Fit(f1, "same")
lat = TLatex()
lat.SetTextSize(0.06)
lat.SetNDC()
lat.SetTextFont(42)
lat.SetTextColor(kBlack)
lat.DrawLatexNDC(0.2, 0.88, f'{name_reso} ({trigger})')
lat.SetTextColor(kRed)
lat.SetTextSize(0.04)
lat.DrawLatexNDC(0.2, 0.82, f'f = p0 + p1/#sqrt{{#it{{N}}_{{trials}}}}')
lat.DrawLatexNDC(0.2, 0.76, f'p0 = {f1.GetParameter(0):.3f} #pm {f1.GetParError(0):.3f}, p1 = {f1.GetParameter(1):.3f} #pm {f1.GetParError(1):.3f}')
canvas.Update()
canvas.SaveAs(f"rysyst_statisticalcomponent_{name_reso}_{trigger}{outlabel}.pdf")
outfile = TFile(f"rysyst_statisticalcomponent_{name_reso}_{trigger}{outlabel}.root", "RECREATE")
for h in htrials:
    h.Write()
gRMS.Write()
f1.Write()
outfile.Close()
input()