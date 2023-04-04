'''
python script to compute inclusive BR(Hb -> Hc + X) reading partial BR from yaml file bHadDecays.yml
'''

import numpy as np
import yaml

with open('bHadDecays.yml', 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

cHadrons = ['D0', 'Dbar0', 'D0D0bar',
            'D+', 'D-', 'D+D-',
            'Ds+', 'Ds-', 'Ds+Ds-',
            'Lc+', 'Lc-', 'Lc+Lc-',
            'D*(2010)+', 'D*(2010)-', 'D*(2010)+D*(2010)-']
BR, BRunc = {}, {}

for cHad in cHadrons:
    for bHad in inputCfg:
        BR[f'{bHad}To{cHad}'] = 0.
        BRunc[f'{bHad}To{cHad}'] = 0.
        for decay in inputCfg[bHad]:
            if cHad in decay:
                BR[f'{bHad}To{cHad}'] += inputCfg[bHad][decay][0]
                BRunc[f'{bHad}To{cHad}'] += inputCfg[bHad][decay][1]**2 # assuming all uncorrelated
                # avoid possible double counting in case of inclusive BRs
                if decay == f'{cHad}X':
                    break
        BRunc[f'{bHad}To{cHad}'] = np.sqrt(BRunc[f'{bHad}To{cHad}'])

for bHad in inputCfg:
    print(f'\n{bHad}')
    for cHad in cHadrons:
        print(f'     {bHad:3} -> {cHad:9} = {BR[f"{bHad}To{cHad}"]:0.4f} +/- {BRunc[f"{bHad}To{cHad}"]:0.4f}')
    print('\n')
