'''
Script with miscellanea utils methods for the analysis
'''

import numpy as np
from ROOT import TH1F, TMath # pylint: disable=import-error,no-name-in-module

def ComputeEfficiency(recoCounts, genCounts, recoCountsError, genCountsError):
    '''
    method to compute efficiency

    Parameters
    ----------
    - recoCounts: number of reconstructed D
    - genCounts: number of genertated D
    - recoCountsError: error on number of reconstructed D
    - genCountsError: error on number of generated D

    Returns
    ----------
    - efficiency, error on efficiency
    '''

    hTmpNum = TH1F('hTmpNum', '', 1, 0, 1)
    hTmpDen = TH1F('hTmpDen', '', 1, 0, 1)
    hTmpNum.SetBinContent(1, recoCounts)
    hTmpDen.SetBinContent(1, genCounts)
    hTmpNum.SetBinError(1, recoCountsError)
    hTmpDen.SetBinError(1, genCountsError)
    hTmpNum.Divide(hTmpNum, hTmpDen, 1., 1, 'B')

    return hTmpNum.GetBinContent(1), hTmpNum.GetBinError(1)


def GetPromptFDYieldsAnalyticMinimisation(effPromptList, effFDList, rawYieldList, \
    effPromptUncList, effFDUncList, rawYieldUncList, precision=1.e-8, nMaxIter=100):
    '''
    method for retrieve prompt and FD corrected yields with an analytic system minimisation

    Parameters
    ----------

    - effPromptList: list of efficiencies for prompt D
    - effFDList: list of efficiencies for FD D
    - rawYieldList: list of raw yields
    - effPromptUncList: list of uncertainties on efficiencies for prompt D
    - effFDUncList: list of uncertainties on efficiencies for FD D
    - rawYieldUncList: list of uncertainties on raw yields
    - precision (float, optional): target precision for minimisation procedure
    - nMaxIter (int, optional): max number of iterations for minimisation procedure

    Returns
    ----------

    - mCorrYield (numpy matrix): corrected yields (Nprompt, NFD)
    - mCovariance (numpy matrix): covariance matrix for corrected yields
    '''

    nCutSets = len(effPromptList)

    mRawYield = np.zeros(shape=(nCutSets, 1))
    mEff = np.zeros(shape=(nCutSets, 2))
    mWeights = np.zeros(shape=(nCutSets, nCutSets))

    mCorrYield = np.zeros(shape=(2, 1))
    mCorrYieldOld = np.zeros(shape=(2, 1))
    mCovariance = np.zeros(shape=(2, 2))

    for iCutSet, (rawYield, effPrompt, effFD) in enumerate(zip(rawYieldList, effPromptList, effFDList)):
        mRawYield.itemset(iCutSet, rawYield)
        mEff.itemset((iCutSet, 0), effPrompt)
        mEff.itemset((iCutSet, 1), effFD)

    mRawYield = np.matrix(mRawYield)
    mEff = np.matrix(mEff)

    for iIter in range(nMaxIter):
        #covariances not taken into account in weight matrix at the moment
        if iIter == 0:
            mCorrYield.itemset(0, 10)
            mCorrYield.itemset(1, 10000)
        for iCutSet, (rawYieldUnc, effPromptUnc, effFDUnc) in enumerate(zip(\
            rawYieldUncList, effPromptUncList, effFDUncList)):
            mWeights.itemset((iCutSet, iCutSet), \
                1. / (rawYieldUnc**2 + effPromptUnc**2 * mCorrYield.item(0) + effFDUnc**2 * mCorrYield.item(1)))

        mWeights = np.matrix(mWeights)
        mEffT = np.transpose(mEff)

        mCovariance = (mEffT * mWeights) * mEff
        mCovariance = np.linalg.inv(mCovariance)
        if mCovariance.item(0, 0) < 1.e-36 or mCovariance.item(1, 1) < 1.e-36:
            print("ERROR: Close to zero or negative diagonal element in covariance matrix:"
                  "likely inversion failed, cannot proceed!")
            return None, None

        mCorrYield = mCovariance * (mEffT * mWeights) * mRawYield

        if (mCorrYield.item(0)-mCorrYieldOld.item(0)) / mCorrYield.item(0) < precision and \
            (mCorrYield.item(1)-mCorrYieldOld.item(1)) / mCorrYield.item(1) < precision:
            break

        mCorrYieldOld = np.copy(mCorrYield)

    return mCorrYield, mCovariance


def SingleGaus(x, par):
    """
    Gaussian function

    Parameters
    ----------

    - x: function variable
    - par: function parameters
        par[0]: normalisation
        par[1]: mean
        par[2]: sigma
    """
    return par[0]*TMath.Gaus(x[0], par[1], par[2], True)


def DoubleGaus(x, par):
    """
    Sum of two Gaussian functions with same mean and different sigma

    Parameters
    ----------

    - x: function variable
    - par: function parameters
        par[0]: normalisation
        par[1]: mean
        par[2]: first sigma
        par[3]: second sigma
        par[4]: fraction of integral in second Gaussian
    """
    firstGaus = TMath.Gaus(x[0], par[1], par[2], True)
    secondGaus = TMath.Gaus(x[0], par[1], par[3], True)
    return par[0] * ((1-par[4])*firstGaus + par[4]*secondGaus)


def DoublePeakSingleGaus(x, par):
    """
    Sum of two Gaussian functions with different mean and sigma

    Parameters
    ----------

    - x: function variable
    - par: function parameters
        par[0]: normalisation first peak
        par[1]: mean first peak
        par[2]: sigma first peak
        par[3]: normalisation second peak
        par[4]: mean second peak
        par[5]: sigma second peak
    """
    firstGaus = par[0]*TMath.Gaus(x[0], par[1], par[2], True)
    secondGaus = par[3]*TMath.Gaus(x[0], par[4], par[5], True)
    return firstGaus + secondGaus


def DoublePeakDoubleGaus(x, par):
    """
    Sum of a double Gaussian function and a single Gaussian function

    Parameters
    ----------

    - x: function variable
    - par: function parameters
        par[0]: normalisation first peak
        par[1]: mean first peak
        par[2]: first sigma first peak
        par[3]: second sigma first peak
        par[4]: fraction of integral in second Gaussian first peak
        par[5]: normalisation second peak
        par[6]: mean second peak
        par[7]: sigma second peak
    """
    firstGaus = TMath.Gaus(x[0], par[1], par[2], True)
    secondGaus = TMath.Gaus(x[0], par[1], par[3], True)
    thirdGaus = par[5]*TMath.Gaus(x[0], par[6], par[7], True)
    return par[0] * ((1-par[4])*firstGaus + par[4]*secondGaus) + thirdGaus
    