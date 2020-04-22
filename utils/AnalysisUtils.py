'''
Script with miscellanea utils methods for the analysis
'''

import ctypes
import numpy as np
import pandas as pd
from ROOT import TH1F, TF1, TMath, TList, TGraphAsymmErrors # pylint: disable=import-error,no-name-in-module

def ComputeEfficiency(recoCounts, genCounts, recoCountsError, genCountsError):
    '''
    Method to compute efficiency

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


# pylint: disable=too-many-locals
def GetPromptFDYieldsAnalyticMinimisation(effPromptList, effFDList, rawYieldList, effPromptUncList, effFDUncList,
                                          rawYieldUncList, corr=True, precision=1.e-8, nMaxIter=100):
    '''
    Method to retrieve prompt and FD corrected yields with an analytic system minimisation

    Parameters
    ----------

    - effPromptList: list of efficiencies for prompt D
    - effFDList: list of efficiencies for FD D
    - rawYieldList: list of raw yields
    - effPromptUncList: list of uncertainties on efficiencies for prompt D
    - effFDUncList: list of uncertainties on efficiencies for FD D
    - rawYieldUncList: list of uncertainties on raw yields
    - corr (bool, optional): whether to compute the correlation
    - precision (float, optional): target precision for minimisation procedure
    - nMaxIter (int, optional): max number of iterations for minimisation procedure

    Returns
    ----------

    - mCorrYield (numpy matrix): corrected yields (Nprompt, NFD)
    - mCovariance (numpy matrix): covariance matrix for corrected yields
    - redChiSquare (float): reduced chi square
    - dicOfMatrices (dictionary): dictionary with all matrices used in minimisation procedure
    '''

    nCutSets = len(effPromptList)

    mRawYield = np.zeros(shape=(nCutSets, 1))
    mEff = np.zeros(shape=(nCutSets, 2))
    mCovSets = np.zeros(shape=(nCutSets, nCutSets))
    mCorrSets = np.zeros(shape=(nCutSets, nCutSets))
    mWeights = np.zeros(shape=(nCutSets, nCutSets))

    mCorrYield = np.zeros(shape=(2, 1))
    mCorrYieldOld = np.zeros(shape=(2, 1))
    mCovariance = np.zeros(shape=(2, 2))
    mRes = np.zeros(shape=(nCutSets, 1))

    for iCutSet, (rawYield, effPrompt, effFD) in enumerate(zip(rawYieldList, effPromptList, effFDList)):
        mRawYield.itemset(iCutSet, rawYield)
        mEff.itemset((iCutSet, 0), effPrompt)
        mEff.itemset((iCutSet, 1), effFD)

    mRawYield = np.matrix(mRawYield)
    mEff = np.matrix(mEff)

    for iIter in range(nMaxIter):
        if iIter == 0:
            mCorrYield.itemset(0, 0)
            mCorrYield.itemset(1, 0)
        for iCutSetRow, (rawYieldUncRow, effPromptUncRow, effFDUncRow) in enumerate(\
            zip(rawYieldUncList, effPromptUncList, effFDUncList)):
            for iCutSetCol, (rawYieldUncCol, effPromptUncCol, effFDUncCol) in enumerate(\
                zip(rawYieldUncList, effPromptUncList, effFDUncList)):
                uncRow = np.sqrt(rawYieldUncRow**2 + effPromptUncRow**2 *
                                 mCorrYield.item(0) + effFDUncRow**2 * mCorrYield.item(1))
                uncCol = np.sqrt(rawYieldUncCol**2 + effPromptUncCol**2 *
                                 mCorrYield.item(0) + effFDUncCol**2 * mCorrYield.item(1))
                if corr and uncRow > 0 and uncCol > 0:
                    if uncRow < uncCol:
                        rho = uncRow / uncCol
                    else:
                        rho = uncCol / uncRow
                else:
                    if iCutSetRow == iCutSetCol:
                        rho = 1.
                    else:
                        rho = 0.
                covRowCol = rho * uncRow * uncCol
                mCovSets.itemset((iCutSetRow, iCutSetCol), covRowCol)
                mCorrSets.itemset((iCutSetRow, iCutSetCol), rho)

        mCovSets = np.matrix(mCovSets)
        mWeights = np.linalg.inv(np.linalg.cholesky(mCovSets))
        mWeights = mWeights.T * mWeights
        mEffT = mEff.T

        mCovariance = (mEffT * mWeights) * mEff
        mCovariance = np.linalg.inv(np.linalg.cholesky(mCovariance))
        mCovariance = mCovariance.T * mCovariance

        mCorrYield = mCovariance * (mEffT * mWeights) * mRawYield
        mRes = mEff * mCorrYield - mRawYield
        mResT = np.transpose(mRes)

        if (mCorrYield.item(0)-mCorrYieldOld.item(0)) / mCorrYield.item(0) < precision and \
            (mCorrYield.item(1)-mCorrYieldOld.item(1)) / mCorrYield.item(1) < precision:
            break

        mCorrYieldOld = np.copy(mCorrYield)

    #reduced chi2
    redChiSquare = mResT * mWeights * mRes / (nCutSets - 2)
    #dictionary with matrices used in minimisation procedure
    dicOfMatrices = {'covMatrix':mCovSets, 'weightMatrix':mWeights, 'corrMatrix':mCorrSets}

    return mCorrYield, mCovariance, float(redChiSquare), dicOfMatrices


def GetPromptFDFractionFc(accEffPrompt, accEffFD, crossSecPrompt, crossSecFD, raaPrompt=1., raaFD=1.):
    '''
    Method to get fraction of prompt / FD fraction with fc method

    Parameters
    ----------

    - accEffPrompt: efficiency times acceptance of prompt D
    - accEffFD: efficiency times acceptance of feed-down D
    - crossSecPrompt: list of production cross sections (cent, min, max) of prompt D in pp collisions from theory
    - crossSecFD: list of production cross sections (cent, min, max) of feed-down D in pp collisions from theory
    - raaPrompt: list of nuclear modification factors (cent, min, max) of prompt D from theory
    - raaFD: list of nuclear modification factors of (cent, min, max) feed-down D from theory

    Returns
    ----------

    - fracPrompt: list of fraction of prompt D (cent, min, max)
    - fracFD: list of fraction of feed-down D (cent, min, max)
    '''
    if not isinstance(crossSecPrompt, list) and isinstance(crossSecPrompt, float):
        crossSecPrompt = [crossSecPrompt]
    if not isinstance(crossSecFD, list) and isinstance(crossSecFD, float):
        crossSecFD = [crossSecFD]
    if not isinstance(raaPrompt, list) and isinstance(raaPrompt, float):
        raaPrompt = [raaPrompt]
    if not isinstance(raaFD, list) and isinstance(raaFD, float):
        raaFD = [raaFD]

    fracPrompt, fracFD = [], []
    if accEffPrompt == 0:
        fracFDCent = 1.
        fracPromptCent = 0.
        fracPrompt = [fracPromptCent, fracPromptCent, fracPromptCent]
        fracFD = [fracFDCent, fracFDCent, fracFDCent]
        return fracPrompt, fracFD
    elif accEffFD == 0:
        fracFDCent = 0.
        fracPromptCent = 1.
        fracPrompt = [fracPromptCent, fracPromptCent, fracPromptCent]
        fracFD = [fracFDCent, fracFDCent, fracFDCent]
        return fracPrompt, fracFD
    else:
        for iSigma, (sigmaP, sigmaF) in enumerate(zip(crossSecPrompt, crossSecFD)):
            for iRaa, (raaP, raaF) in enumerate(zip(raaPrompt, raaFD)):
                if iSigma == 0 and iRaa == 0:
                    fracPromptCent = 1./(1 + accEffFD / accEffPrompt * sigmaF / sigmaP * raaF / raaP)
                    fracFDCent = 1./(1 + accEffPrompt / accEffFD * sigmaP / sigmaF * raaP / raaF)
                else:
                    fracPrompt.append(1./(1 + accEffFD / accEffPrompt * sigmaF / sigmaP * raaF / raaP))
                    fracFD.append(1./(1 + accEffPrompt / accEffFD * sigmaP / sigmaF * raaP / raaF))

    if fracPrompt and fracFD:
        fracPrompt.sort()
        fracFD.sort()
        fracPrompt = [fracPromptCent, fracPrompt[0], fracPrompt[-1]]
        fracFD = [fracFDCent, fracFD[0], fracFD[-1]]
    else:
        fracPrompt = [fracPromptCent, fracPromptCent, fracPromptCent]
        fracFD = [fracFDCent, fracFDCent, fracFDCent]

    return fracPrompt, fracFD


# pylint: disable=too-many-arguments, too-many-branches
def GetFractionNb(rawYield, accEffSame, accEffOther, crossSec, deltaPt, deltaY, BR, nEvents, \
    sigmaMB, raaRatio=1., taa=1., ppRef=1.):
    '''
    Method to get fraction of prompt / FD fraction with Nb method

    Parameters
    ----------

    - accEffSame: efficiency times acceptance of prompt (feed-down) D
    - accEffOther: efficiency times acceptance of feed-down (prompt) D
    - crossSec: list of production cross sections (cent, min, max) of feed-down (prompt)
      D in pp collisions from theory
    - deltaPt: width of pT interval
    - deltaY: width of Y interval
    - BR: branching ratio for the chosen decay channel
    - nEvents: number of events corresponding to the raw yields
    - sigmaMB: MB cross section (=1 for p-Pb and Pb-Pb)
    - raaRatio: list of D nuclear modification factor ratios
      feed-down / prompt (prompt / feed-down) (cent, min, max) (=1 in case of pp)
    - taa: average nuclear overlap function (=1 in case of pp)
    - ppRef: value of pp reference for prompt (feed-down) D (=1 in case of pp)

    Returns
    ----------

    - frac: list of fraction of prompt (feed-down) D (cent, min, max)
    '''
    if not isinstance(crossSec, list) and isinstance(crossSec, float):
        crossSec = [crossSec]

    if not isinstance(raaRatio, list) and isinstance(raaRatio, float):
        raaRatio = [raaRatio]

    frac = []
    for iSigma, sigma in enumerate(crossSec):
        for iRaaRatio, raaRat in enumerate(raaRatio):
            raaOther = 1.
            if iSigma == 0 and iRaaRatio == 0:
                if raaRat == 1. and ppRef == 1. and taa == 1.: #pp
                    fracCent = 1 - sigma * deltaPt * deltaY * accEffOther * BR * nEvents * 2 / rawYield / sigmaMB
                else: #p-Pb or Pb-Pb: iterative evaluation of Raa needed
                    deltaRaa = 1.
                    while deltaRaa > 1.e-3:
                        fracCent = 1 - taa * raaRat * raaOther * sigma * \
                            deltaPt * deltaY * accEffOther * BR * nEvents * 2 / rawYield
                        raaOtherOld = raaOther
                        raaOther = fracCent * rawYield * sigmaMB / (2 * accEffSame * deltaPt * deltaY * BR * nEvents)
                        deltaRaa = abs((raaOther-raaOtherOld) / raaOther)

            else:
                if raaRat == 1. and ppRef == 1. and taa == 1.: #pp
                    frac.append(1 - sigma * deltaPt * deltaY * accEffOther * BR * nEvents * 2 / rawYield / sigmaMB)
                else:
                    deltaRaa = 1.
                    fracTmp = 1.
                    while deltaRaa > 1.e-3:
                        fracTmp = 1 - taa * raaRat * raaOther * sigma * \
                            deltaPt * deltaY * accEffOther * BR * nEvents * 2 / rawYield
                        raaOtherOld = raaOther
                        raaOther = fracTmp * rawYield * sigmaMB / (2 * accEffSame * deltaPt * deltaY * BR * nEvents)
                        deltaRaa = abs((raaOther-raaOtherOld) / raaOther)
                    frac.append(fracTmp)

    if frac:
        frac.sort()
        frac = [fracCent, frac[0], frac[-1]]
    else:
        frac = [fracCent, fracCent, fracCent]

    return frac


def SingleGaus(x, par):
    '''
    Gaussian function

    Parameters
    ----------

    - x: function variable
    - par: function parameters
        par[0]: normalisation
        par[1]: mean
        par[2]: sigma
    '''
    return par[0]*TMath.Gaus(x[0], par[1], par[2], True)


def DoubleGaus(x, par):
    '''
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
    '''
    firstGaus = TMath.Gaus(x[0], par[1], par[2], True)
    secondGaus = TMath.Gaus(x[0], par[1], par[3], True)
    return par[0] * ((1-par[4])*firstGaus + par[4]*secondGaus)


def DoublePeakSingleGaus(x, par):
    '''
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
    '''
    firstGaus = par[0]*TMath.Gaus(x[0], par[1], par[2], True)
    secondGaus = par[3]*TMath.Gaus(x[0], par[4], par[5], True)
    return firstGaus + secondGaus


def DoublePeakDoubleGaus(x, par):
    '''
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
    '''
    firstGaus = TMath.Gaus(x[0], par[1], par[2], True)
    secondGaus = TMath.Gaus(x[0], par[1], par[3], True)
    thirdGaus = par[5]*TMath.Gaus(x[0], par[6], par[7], True)
    return par[0] * ((1-par[4])*firstGaus + par[4]*secondGaus) + thirdGaus


def GetExpectedBkgFromSideBands(hMassData, bkgFunc='pol2', nSigmaForSB=4, hMassSignal=None, mean=-1., sigma=-1.,
                                hMassSecPeak=None, meanSecPeak=-1., sigmaSecPeak=-1.):
    '''
    Helper method to get the expected bkg from side-bands

    Parameters
    ----------

    - hMassData: invariant-mass histogram from which extract the estimated bkg
    - bkgFunc: expression for bkg fit function
    - nSigmaForSB: number of sigmas away from the invariant-mass peak to define SB windows
    - hMassSignal: invariant-mass histogram for the signal used to get mean and sigma
                   not needed in case of passed mean and sigma parameters
    - mean: mean of invariant-mass peak of the signal
                   not needed in case of passed hMassSignal
    - sigma: width of invariant-mass peak of the signal
                   not needed in case of passed hMassSignal
    - hMassSecPeak: invariant-mass histogram for the second peak (only Ds) used to get meanSecPeak and sigmaSecPeak
                   not needed in case of passed meanSecPeak and sigmaSecPeak parameters
    - meanSecPeak: mean of invariant-mass peak of the second peak (only Ds)
                   not needed in case of passed hMassSecPeak
    - sigmaSecPeak: width of invariant-mass peak of the second peak (only Ds)
                   not needed in case of passed hMassSecPeak

    Returns
    ----------

    - expBkg3s: expected background within 3 sigma from signal peak mean
    - hMassData: SB histogram with fit function
    '''
    if hMassSignal:
        funcSignal = TF1('funcSignal', SingleGaus, 1.6, 2.2, 3)
        funcSignal.SetParameters(
            hMassSignal.Integral() * hMassSignal.GetBinWidth(1), hMassSignal.GetMean(), hMassSignal.GetRMS())
        hMassSignal.Fit('funcSignal', 'Q0')
        mean = funcSignal.GetParameter(1)
        sigma = funcSignal.GetParameter(2)
    if hMassSecPeak:
        funcSignal.SetParameters(
            hMassSecPeak.Integral() * hMassSecPeak.GetBinWidth(1), hMassSecPeak.GetMean(), hMassSecPeak.GetRMS())
        hMassSecPeak.Fit('funcSignal', 'Q0')
        meanSecPeak = funcSignal.GetParameter(1)
        sigmaSecPeak = funcSignal.GetParameter(2)

    for iMassBin in range(1, hMassData.GetNbinsX()+1):
        massLowLimit = hMassData.GetBinLowEdge(iMassBin)
        massUpLimit = hMassData.GetBinLowEdge(iMassBin) + hMassData.GetBinWidth(iMassBin)

        if massLowLimit > mean - nSigmaForSB * sigma and massUpLimit < mean + nSigmaForSB * sigma:
            hMassData.SetBinContent(iMassBin, 0.)
            hMassData.SetBinError(iMassBin, 0.)
        elif meanSecPeak > 0 and sigmaSecPeak > 0:
            if massLowLimit > meanSecPeak - nSigmaForSB * sigmaSecPeak and \
                massUpLimit < meanSecPeak + nSigmaForSB * sigmaSecPeak:
                hMassData.SetBinContent(iMassBin, 0.)
                hMassData.SetBinError(iMassBin, 0.)

    funcBkg = TF1('funcBkg', bkgFunc, 1.6, 2.2)
    hMassData.Fit(funcBkg, 'Q')
    expBkg3s = funcBkg.Integral(mean - 3 * sigma, mean + 3 * sigma) / hMassData.GetBinWidth(1)
    return expBkg3s, hMassData


def GetExpectedBkgFromMC(hMassBkg, hMassSignal=None, mean=-1., sigma=-1., doFit=True, bkgFunc='pol3'):
    '''
    Helper method to get the expected bkg from MC

    Parameters
    ----------

    - hMassBkg: invariant-mass histogram of background
    - hMassSignal: invariant-mass histogram for the signal used to get mean and sigma
                   not needed in case of passed mean and sigma parameters
    - mean: mean of invariant-mass peak of the signal
                   not needed in case of passed hMassSignal
    - sigma: width of invariant-mass peak of the signal
                   not needed in case of passed hMassSignal
    - doFit: flag to enable fit of bkg distribution (useful with low stat)
    - bkgFunc: expression for bkg fit function (if fit enabled)

    Returns
    ----------

    - expBkg3s: expected background within 3 sigma from signal peak mean
    '''
    if hMassSignal:
        funcSignal = TF1('funcSignal', SingleGaus, 1.6, 2.2, 3)
        funcSignal.SetParameters(hMassSignal.Integral() * hMassSignal.GetBinWidth(1),
                                 hMassSignal.GetMean(), hMassSignal.GetRMS())
        hMassSignal.Fit('funcSignal', 'Q0')
        mean = funcSignal.GetParameter(1)
        sigma = funcSignal.GetParameter(2)

    if doFit:
        funcBkg = TF1('funcBkg', bkgFunc, 1.6, 2.2)
        hMassBkg.Fit(funcBkg, 'Q')
        hMassBkg.Write()
        expBkg3s = funcBkg.Integral(mean - 3 * sigma, mean + 3 * sigma) / hMassBkg.GetBinWidth(1)
    else:
        massBinMin = hMassBkg.GetXaxis().FindBin(mean - 3 * sigma)
        massBinMax = hMassBkg.GetXaxis().FindBin(mean + 3 * sigma)
        expBkg3s = hMassBkg.Integral(massBinMin, massBinMax)

    return expBkg3s


def GetExpectedSignal(crossSec, deltaPt, deltaY, effTimesAcc, frac, BR, fractoD, nEv, sigmaMB=1, TAA=1, RAA=1):
    '''
    Helper method to get expected signal from MC and predictions

    Parameters
    ----------

    - crossSec: prediction for differential cross section in pp
    - deltaPt: pT interval
    - deltaY: Y interval
    - effTimesAcc: efficiency times acceptance for prompt or feed-down
    - frac: either prompt or feed-down fraction
    - BR: branching ratio of the decay channel
    - fracToD: fragmentation fraction
    - nEv: number of expected events
    - sigmaMB: hadronic cross section for MB
    - TAA: average overlap nuclear function
    - RAA: expected nuclear modification factor

    Returns
    ----------

    - expected signal
    '''

    return 2 * crossSec * deltaPt * deltaY * effTimesAcc * BR * fractoD * nEv * TAA * RAA / frac / sigmaMB


def ComputeCrossSection(rawY, uncRawY, frac, uncFrac, effTimesAcc, deltaPt, deltaY, sigmaMB, nEv, BR):
    '''
    Method to compute cross section and its statistical uncertainty
    Only the statistical uncertainty on the raw yield and prompt (feed-down)
    fraction are considered (the others are systematics)

    Parameters
    ----------

    - rawY: raw yield
    - uncRawY: raw-yield statistical uncertainty
    - frac: either prompt or feed-down fraction
    - uncFrac: uncertainty on prompt or feed-down fraction
    - effTimesAcc: efficiency times acceptance for prompt or feed-down
    - deltaPt: pT interval
    - deltaY: Y interval
    - sigmaMB: hadronic cross section for MB
    - nEv: number of events
    - BR: branching ratio of the decay channel


    Returns
    ----------

    - crossSection: cross section
    - crossSecUnc: cross-section statistical uncertainty
    '''

    crossSection = rawY * frac * sigmaMB / (2 * deltaPt * deltaY * effTimesAcc * nEv * BR)
    crossSecUnc = np.sqrt((uncRawY / rawY)**2 + (uncFrac / frac)**2) * crossSection

    return crossSection, crossSecUnc


def MergeHists(listOfHists):
    '''
    Method to merge histos

    Parameters
    ----------

    - listOfHists: python list of histos

    Returns
    ----------

    - hMerged: merged histo

    '''
    listMerge = TList()
    for iHist, hist in enumerate(listOfHists):
        if iHist == 0:
            hMerged = hist.Clone()
        else:
            listMerge.Add(hist)
    hMerged.Merge(listMerge)
    return hMerged


def ApplySplineFuncToColumn(df, column, spline, minRange=-1.e10, maxRange=1.e10):
    '''
    Method to apply a function to a pandas column via a spline object

    Parameters
    ----------

    - df: input pandas.Dataframe
    - column: column of the pandas dataframe to which apply the spline
    - spline: spline (scipy.interpolate.InterpolatedUnivariateSpline object)
    - minRange: minimum of the range of the user-defined validity of the spline
    - maxRange: maximum of the range of the user-defined validity of the spline

    Returns
    ----------

    - y: pandas.Series with result of the function application to column

    '''

    y = []
    for x in df[column].values:
        if minRange <= x <= maxRange:
            y.append(spline(x))
        elif x < minRange:
            y.append(spline(minRange))
        else:
            y.append(spline(maxRange))

    y = pd.Series(y)

    return y


def ComputeRatioDiffBins(hNum, hDen, uncOpt=''):
    '''
    Method to compute ratio between histograms with different bins (but compatible)

    Parameters
    ----------

    - hNum: histogram for numerator
    - hDen: histogram for denominator
    - uncOpt: uncertainty option as in ROOT.TH1.Divide

    Returns
    ----------

    - hRatio: ratio histogram

    '''

    ptMinNum = hNum.GetBinLowEdge(1)
    ptMaxNum = hNum.GetXaxis().GetBinUpEdge(hNum.GetNbinsX())
    ptMinDen = hDen.GetBinLowEdge(1)
    ptMaxDen = hDen.GetXaxis().GetBinUpEdge(hDen.GetNbinsX())
    if ptMinNum < ptMinDen:
        ptMin = ptMinDen
    else:
        ptMin = ptMinNum
    if ptMaxNum > ptMaxDen:
        ptMax = ptMaxDen
    else:
        ptMax = ptMaxNum

    if hNum.GetNbinsX() < hDen.GetNbinsX():
        ptLimsRatio = np.array(hNum.GetXaxis().GetXbins(), 'd')
    else:
        ptLimsRatio = np.array(hDen.GetXaxis().GetXbins(), 'd')
    ptLimsRatio = ptLimsRatio[(ptLimsRatio >= ptMin) & (ptLimsRatio <= ptMax)]
    nPtBins = len(ptLimsRatio)-1

    hRatio = TH1F('hRatio', f';{hNum.GetXaxis().GetTitle()};ratio', nPtBins, ptLimsRatio)
    hNumReb = TH1F('hNumReb', '', nPtBins, ptLimsRatio)
    hDenReb = TH1F('hDenReb', '', nPtBins, ptLimsRatio)

    for iPtRatio in range(1, hRatio.GetNbinsX()+1):
        deltaPt = ptLimsRatio[iPtRatio]-ptLimsRatio[iPtRatio-1]
        num, numUnc, den, denUnc = (0 for _ in range(4))
        for iPtNum in range(1, hNum.GetNbinsX()+1):
            if hNum.GetBinLowEdge(iPtNum) >= ptLimsRatio[iPtRatio-1] and \
                hNum.GetXaxis().GetBinUpEdge(iPtNum) <= ptLimsRatio[iPtRatio]:
                num += hNum.GetBinContent(iPtNum) * hNum.GetBinWidth(iPtNum)
                numUnc += hNum.GetBinError(iPtNum)**2 * hNum.GetBinWidth(iPtNum)**2 # considered uncorrelated
        hNumReb.SetBinContent(iPtRatio, num/deltaPt)
        hNumReb.SetBinError(iPtRatio, np.sqrt(numUnc)/deltaPt)
        for iPtDen in range(1, hDen.GetNbinsX()+1):
            if hDen.GetBinLowEdge(iPtDen) >= ptLimsRatio[iPtRatio-1] and \
                hDen.GetXaxis().GetBinUpEdge(iPtDen) <= ptLimsRatio[iPtRatio]:
                den += hDen.GetBinContent(iPtDen) * hDen.GetBinWidth(iPtDen)
                denUnc += hDen.GetBinError(iPtDen)**2 * hDen.GetBinWidth(iPtDen)**2 # considered uncorrelated
        hDenReb.SetBinContent(iPtRatio, den/deltaPt)
        hDenReb.SetBinError(iPtRatio, np.sqrt(denUnc)/deltaPt)

    hRatio.Divide(hNumReb, hDenReb, 1., 1., uncOpt)

    return hRatio


def ScaleGraph(graph, scaleFactor):
    '''
    Helper method to scale a TGraph

    Parameters
    ----------

    - graph: graph to scale
    - scaleFactor: scale factor
    '''
    for iPt in range(graph.GetN()):
        x, y = ctypes.c_double(), ctypes.c_double()
        graph.GetPoint(iPt, x, y)
        graph.SetPoint(iPt, x.value, y.value * scaleFactor)
        yUncLow = graph.GetErrorYlow(iPt)
        yUncHigh = graph.GetErrorYhigh(iPt)
        graph.SetPointEYlow(iPt, yUncLow * scaleFactor)
        graph.SetPointEYhigh(iPt, yUncHigh * scaleFactor)


def DivideGraphByHisto(gNum, hDen, useHistoUnc=True):
    '''
    Helper method to divide a TGraph by a TH1 (assuming same binning)

    Parameters
    ----------

    - gNum: graph to divide (numerator)
    - hDen: histogram (denominator)

    Returns
    ----------

    - gRatio: resulting graph
    '''
    if gNum.GetN() != hDen.GetNbinsX():
        print('ERROR: only graphs and histos with same number of bins can be divided!')
        return None

    gRatio = TGraphAsymmErrors(0)
    for iPt in range(gNum.GetN()):
        x, num = ctypes.c_double(), ctypes.c_double()
        gNum.GetPoint(iPt, x, num)
        xUncLow = gNum.GetErrorXlow(iPt)
        xUncHigh = gNum.GetErrorXhigh(iPt)
        numUncLow = gNum.GetErrorYlow(iPt)
        numUncHigh = gNum.GetErrorYhigh(iPt)
        den = hDen.GetBinContent(iPt+1)
        if useHistoUnc:
            ratioUncLow = np.sqrt((numUncLow/num)**2 + (hDen.GetBinError(iPt+1)/den)**2) * num/den
            ratioUncHigh = np.sqrt((numUncHigh/num)**2 + (hDen.GetBinError(iPt+1)/den)**2) * num/den
        else:
            ratioUncLow = numUncLow/num.value * num.value/den
            ratioUncHigh = numUncHigh/num.value * num.value/den
        gRatio.SetPoint(iPt, x.value, num.value/den)
        gRatio.SetPointError(iPt, xUncLow, xUncHigh, ratioUncLow, ratioUncHigh)

    return gRatio

def ApplyVariationToList(listToVary, relVar, option='decreasing'):  
    '''
    Helper method to apply a relative variation to a list of numbers

    Parameters
    ----------

    - listToVary: list of values to be varied
    - relVar: relative variation
    - option: option for variation among
        - upshift: all the values are shifted upwards by relVar
        - downshift: all the values are shifted downwards by relVar
        - decreasing: the first value is not varied; the next ones are decreased smoothly up to relVar for the last one
        - increasing: the first value is not varied; the next ones are increased smoothly up to relVar for the last one

    Returns
    ----------

    - listVaried: list of varied values
    '''

    if option not in ['upshift', 'downshift', 'decreasing', 'increasing']:
        print(f'ERROR: option for variation of list not valid! Returning None')
        return None

    if option == 'upshift':
        listVaried = [el + el*relVar for el in listToVary]
    elif option == 'downshift':
        listVaried = [el - el*relVar for el in listToVary]
    elif option == 'decreasing':
        listVaried = [el - el*relVar/len(listToVary)*(iEl+1) for iEl, el in enumerate(listToVary)]
    elif option == 'increasing':
        listVaried = [el + el*relVar/len(listToVary)*(iEl+1) for iEl, el in enumerate(listToVary)]

    return listVaried
