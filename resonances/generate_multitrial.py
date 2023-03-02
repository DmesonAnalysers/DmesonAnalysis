'''
python script to generate the multitrial ditribution for systematic checks
run: python generate_multitrial.py
'''
import sys
from ROOT import TFile, TH1F

pdg_reso = 435
trigger = "HM"

if pdg_reso == 10433:
    name_reso = "Ds1plus"
    label_mass = r"M($\mathrm{D}^{*+}\mathrm{K_S^0}$)$-$M($\mathrm{D}^{*+}$)"
    pdg_v0 = 310
    pdg_d = 413
    signal_limits = [0.51, 0.53]
elif pdg_reso == 435:
    name_reso = "Ds2starplus"
    label_mass = r"M($\mathrm{D}^+\mathrm{K_S^0}$)$-$M($\mathrm{D}^+$)"
    pdg_v0 = 310
    pdg_d = 411
    signal_limits = [0.65, 0.75]
else:
    print(f"ERROR: pdg code {pdg_reso} not supported")
    sys.exit()

#______________________________________________________________________________
# Load reference file
ref_file = f'/home/stefano/Desktop/cernbox/Ds_resonances/raw_yield/mass_{name_reso}_pt2.0-24.0_{trigger}.root'
print(f'\033[1m\033[92m Loading input file {ref_file} \033[0m')
infile_ref = TFile(ref_file, "READ")
h_data = infile_ref.Get("hdata")
hbkg_ref = infile_ref.Get("bkg_0") # non-trustable normalisation
hsignal_ref = infile_ref.Get("signal_0") # non-trustable normalisation
hsignal = infile_ref.Get("h_rawyields")
hbkg_all = infile_ref.Get("h_bkg_allrange")
nsignal = hsignal.GetBinContent(1)
nbkg = hbkg_all.GetBinContent(1)

#______________________________________________________________________________
# Generate multi-trial distribution
outfile = TFile(f"multitrialcheck_{name_reso}_{trigger}.root", "RECREATE")
for i in range(100):
    htrial = TH1F()
    htrial = h_data.Clone(f"htrial_{i}")
    htrial.Reset()
    htrial.SetName(f"htrial_{i}")
    htrial.FillRandom(hsignal_ref, int(nsignal))
    htrial.FillRandom(hbkg_ref, int(nbkg))
    htrial.Write()
print(f"Output file {outfile.GetName()} saved")