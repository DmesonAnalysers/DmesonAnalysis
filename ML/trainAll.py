'''
Script to parallelize trainings and application.
run: python trainAll.py configFileName.yml

'''
from multiprocessing import Pool
import sys
import MLClassification
import argparse
import yaml

sys.path.append('..')
from utils.Misc import exec_time

def trainAndApply(traincfgfile, applycfgfile):
    MLClassification.main(traincfgfile, doTraining= True, doApplication=False)
    MLClassification.main(applycfgfile, doTraining= False, doApplication=True)

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileNameML.yml', help='config file name for ml')
    args = parser.parse_args()

    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    trainingList = inputCfg["cfgtrainfiles"]
    applicationList = inputCfg["cfgapplyfiles"]

    with Pool(inputCfg['njobs']) as p:
        print(p.starmap(trainAndApply, tuple(zip(trainingList, applicationList))))

main()