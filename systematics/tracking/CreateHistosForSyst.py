'''
python script for the creation of the input file with tracking systematics per track
run: python CreateHistosForSyst.py period [--Ds] [--Dplus]
'''

import sys
import argparse
import yaml
import numpy as np
from ROOT import TH1F, TFile # pylint: disable=import-error,no-name-in-module

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('period', metavar='text', default='LHC17pq', help='data period for systematic evaluation')
parser.add_argument('--Dplus', action='store_true', default=False, help='enable comparison for D+')
parser.add_argument('--Ds', action='store_true', default=False, help='enable comparison for Ds')
args = parser.parse_args()

if not args.Dplus and not args.Ds:
    print('ERROR: you should enable the syst uncertainty evaluation for either D+ or Ds! Exit')
    sys.exit()
elif args.Dplus and args.Ds:
    print('ERROR: you should enable the syst uncertainty evaluation for either D+ or Ds! Exit')
    sys.exit()
elif args.Dplus:
    meson = 'Dplus'
else:
    meson = 'Ds'

#config with input file details
with open('db_tracking_syst_singletrack.yml', 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

ptMins = inputCfg[args.period]['TPCquality'][meson]['ptmins']
ptMaxs = inputCfg[args.period]['TPCquality'][meson]['ptmaxs']
ptLims = ptMins.copy()
ptLims.append(ptMaxs[-1])
systTPC = inputCfg[args.period]['TPCquality'][meson]['syst']

ptMinsME = inputCfg[args.period]['ITSTPCme']['ptmins']
ptMaxsME = inputCfg[args.period]['ITSTPCme']['ptmaxs']
ptLimsME = ptMinsME.copy()
ptLimsME.append(ptMaxsME[-1])
systME = inputCfg[args.period]['ITSTPCme']['syst']

hTrk = TH1F('hTrEff', ';#it{p}_{T} (GeV/#it{c});rel syst', len(ptMins), np.array(ptLims, 'f'))
for iPt, syst in enumerate(systTPC):
    hTrk.SetBinContent(iPt+1, syst)

hME = TH1F('hME', ';#it{p}_{T} (GeV/#it{c});rel syst', len(ptMinsME), np.array(ptLimsME, 'f'))
for iPt, syst in enumerate(systME):
    hME.SetBinContent(iPt+1, syst)

outfile = TFile(f'singletracksyst/traking_ME_piK_syst_{meson}_{args.period}.root', 'recreate')
hTrk.Write('hTrEff')
hME.Write('h')
outfile.Close()
