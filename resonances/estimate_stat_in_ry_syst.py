'''
python script to estimate the statistical uncertainty in the systematic check
run: python estimate_stat_in_ry_syst.py
'''
import sys
import numpy as np
from ROOT import TFile, TH1F, TGraphErrors, kRed, TF1, kBlack, TCanvas, TLatex, kBlue
sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
SetGlobalStyle(padtopmargin=0.05, padleftmargin=0.18,padbottommargin=0.15,
               palette=55, labelsize=0.04, titlesize=0.06,
               labeloffset=0.008, titleoffsety=1.55, titleoffsetx=1., titleoffsetz=1.,
               opttitle=0, optstat=0)
name_reso = "Ds2starplus"
trigger = "MB"
outlabel = "_FINAL"

#______________________________________________________________________________
# Load reference file
mults = [1, 20, 40, 100] if trigger == "MB" else [1, 10, 20, 40, 100] # at the moment exluding 10 for MB
gRMS = TGraphErrors()
SetObjectStyle(gRMS, color=kBlack, markerstyle=20, markersize=0.5)
rms, rms_unc = [], []
hMult = []
for i, mult in enumerate(mults):
    meanval = 0.
    meanval_unc = 0.
    hMult.append(TH1F(f"hMult{mult}", f"hMult{i}", 50, 0, 0.5))
    hMult[-1].SetDirectory(0)
    for i in range(10):
        if trigger == "HM" and i in [4]: # skip bad trials in HM (4)
            continue
        file = TFile(f'./Multitrial_{name_reso}_pt2.0-24.0_{trigger}_test_bash_htrial{i}_multfact{mult:.1f}.root')
        histo = file.Get('h_syst')
        sigma = histo.GetRMS()
        mean = histo.GetMean()
        sigma_unc = histo.GetRMSError()
        mean_unc = histo.GetMeanError()
        hMult[-1].Fill(sigma/mean)
    rms.append(hMult[-1].GetMean())
    rms_unc.append(hMult[-1].GetRMS()/mult)

canvas = TCanvas("canvas", "canvas", 800, 800)
ymax = max(rms) + 0.1
ymin = min(rms) - 0.02
hFrame = canvas.DrawFrame(0, ymin, max(mults)+5, ymax, ";N_{evt}/N_{evt}(data);Average RMS/Mean")
for i, (mult, rms, rms_unc) in enumerate(zip(mults, rms, rms_unc)):
    gRMS.SetPoint(i, mult, rms)
    gRMS.SetPointError(i, 0, rms_unc)
yhighfit = 30. if trigger == "HM" else 100.

for mult in mults:
    f1 = TF1(f"f1_multfact{mult}", "[0] + [1]/sqrt(x)", 1.e-2, mult+2)
    gRMS.Fit(f1, "same r")
    chi2 = f1.GetChisquare()
    ndf = f1.GetNDF()
    if ndf != 0:
        red_chi2 = chi2/ndf
        if red_chi2 > 5: # keep extending the fit range until the chi2/ndf is reasonable
            print(f"\033[1;31mWARNING: red_chi2 = {red_chi2:.2f} > 2\033[0m")
            continue
        else:
            print(f"red_chi2 = {red_chi2:.2f}")
            yhighfit = mult+2
f1 = TF1("f1", "[0] + [1]/sqrt(x)", 1.e-2, yhighfit)
f1.SetLineColor(kRed+1)
f1.SetLineStyle(2)
gRMS.Draw('PZ same')
gRMS.Fit(f1, "same r")
lat = TLatex()
lat.SetTextSize(0.06)
lat.SetNDC()
lat.SetTextFont(42)
lat.SetTextColor(kBlack)
lat.DrawLatexNDC(0.2, 0.88, f'{name_reso} ({trigger})')
lat.SetTextColor(kRed+1)
lat.SetTextSize(0.04)
lat.DrawLatexNDC(0.2, 0.82, '#it{f} = p0 + p1/#sqrt{N_{evt}}')
lat.DrawLatexNDC(0.2, 0.76, f'p0 = {f1.GetParameter(0):.3f} +- {f1.GetParError(0):.3f}, p1 = {f1.GetParameter(1):.3f} #pm {f1.GetParError(1):.3f}')
lat.DrawLatexNDC(0.2, 0.20, f'#chi^{{2}}/ndf = {f1.GetChisquare()/f1.GetNDF():.2f}')
canvas.Update()
print(f'\033[1;32mINFO: RMS/Mean = {f1.GetParameter(1):.3f} +- {f1.GetParError(1):.3f}\033[0m')
print(f'\033[1;32mINFO: SYST = {f1.GetParameter(0):.3f} +- {f1.GetParError(0):.3f}\033[0m')
canvas.SaveAs(f"rysyst_statisticalcomponent_{name_reso}_{trigger}{outlabel}.pdf")
outfile = TFile(f"rysyst_statisticalcomponent_{name_reso}_{trigger}{outlabel}.root", "RECREATE")
for h in hMult:
    h.Write()
gRMS.Write()
f1.Write()
canvas.Write()
outfile.Close()
input()
