'''
python script to generate the multitrial ditribution for systematic checks
run: python generate_multitrial.py reffile.root multitrialcheck.root mult_factor
'''
import sys
import argparse
from ROOT import TFile, TH1F

parser = argparse.ArgumentParser(description="Arguments")
parser.add_argument("reffile", metavar="text",
                    default="reffile.root", help="input file")
parser.add_argument('outfilename', metavar=('text'),
                    default="multitrialcheck.root", help='Output file name')
parser.add_argument('mult_factor', type=float,
                    default=1., help='multiplicative factor for signal and background')
args = parser.parse_args()

if "Ds1" in args.reffile:
    name_reso = "Ds1plus"
    label_mass = r"M($\mathrm{D}^{*+}\mathrm{K_S^0}$)$-$M($\mathrm{D}^{*+}$)"
    pdg_v0 = 310
    pdg_d = 413
elif "Ds2" in args.reffile:
    name_reso = "Ds2starplus"
    label_mass = r"M($\mathrm{D}^+\mathrm{K_S^0}$)$-$M($\mathrm{D}^+$)"
    pdg_v0 = 310
    pdg_d = 411
else:
    print(f"ERROR: resonance not recognized, only Ds1 and Ds2 are supported")
    sys.exit()

if "MB" in args.reffile:
    trigger = "MB"
elif "HM" in args.reffile:
    trigger = "HM"
else:
    print(f"ERROR: trigger not recognized, only MB and HM are supported")
    sys.exit()

#______________________________________________________________________________
# Load reference file
ref_file = args.reffile
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
outfilename = args.outfilename.replace(".root", f"_multfactor{args.mult_factor}.root")
outfile = TFile(f"{outfilename}", "RECREATE")
for i in range(10):
    htrial = TH1F()
    htrial = h_data.Clone(f"htrial{i}_multfact{args.mult_factor}")
    htrial.Reset()
    htrial.FillRandom(hsignal_ref, int(nsignal*args.mult_factor)) # TODO: move to unbinned case
    htrial.FillRandom(hbkg_ref, int(nbkg*args.mult_factor))
    htrial.Write()
print(f"Output file {outfile.GetName()} saved")