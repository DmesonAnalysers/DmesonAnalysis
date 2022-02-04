import sys
import numpy as np # se importo root devo anche importare numpy (???)
from ROOT import TFile
sys.path.append('..')
from utils.DfUtils import GetObjectFromFile #pylint: disable=wrong-import-position,import-error

fn = '../../CharmingAnalyses/DKDpi/oton_mctruth/data/AnalysisResults_all.root'
objn = 'HM_CharmFemto_Dkaon_Results0/HM_CharmFemto_Dkaon_Results0/Particle0_Particle2/MEDist_Particle0_Particle2'

o = GetObjectFromFile(fn, objn)

print(o)