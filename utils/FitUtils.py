'''
Module with function definitions and fit utils
'''

from ROOT import TMath, TF1, kBlue, kGreen # pylint: disable=import-error,no-name-in-module

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


def VoigtFunc(x, par):
    '''
    Voigtian function

    Parameters
    ----------

    - x: function variable
    - par: function parameters
        par[0]: normalisation
        par[1]: mean
        par[2]: sigma
        par[3]: gamma
    '''

    return par[0] * TMath.Voigt(x[0]-par[1], par[2], par[3])


def ExpoPowLaw(x, par):
    '''
    Exponential times power law function

    Parameters
    ----------

    - x: function variable
    - par: function parameters
        par[0]: normalisation
        par[1]: mass (lowest possible value)
        par[2]: expo slope
    '''

    return par[0] * TMath.Sqrt(x[0] - par[1]) * TMath.Exp(-1. * par[2] * (x[0] - par[1]))


# pylint: disable=too-many-instance-attributes
class BkgFuncCreator:
    '''
    Class to handle custom background functions as done by AliHFInvMassFitter. Mainly designed
    to provide functions for sidebands fitting
    '''

    def __init__(self, funcName, minMassHisto, maxMassHisto, numSigmaSideBands=0., peakMass=0.,
                 peakSigma=0., secPeakMass=0., secPeakSigma=0.):
        self.funcName = funcName
        self.minMass = minMassHisto
        self.maxMass = maxMassHisto
        self.peakMass = peakMass
        self.peakDelta = peakSigma * numSigmaSideBands
        self.secPeakMass = secPeakMass
        self.secPeakDelta = secPeakSigma * numSigmaSideBands

        self.removePeak = False
        self.removeSecPeak = False
        if self.peakMass > 0. and self.peakDelta > 0.:
            self.removePeak = True
        if self.secPeakMass > 0. and self.secPeakDelta > 0.:
            self.removeSecPeak = True

    def _ExpIntegralNorm(self, x, par):
        '''
        Exponential function normalized to its integral.
        See AliHFInvMassFitter::FitFunction4Bkg for more information.

        Parameters
        ----------

        - x: function variable
        - par: function parameters
            par[0]: normalisation (integral of background)
            par[1]: expo slope
        '''
        norm = par[0] * par[1] / (TMath.Exp(par[1] * self.maxMass) - TMath.Exp(par[1] * self.minMass))
        return norm * TMath.Exp(par[1] * x[0])

    def _SideBandsFunc(self, x, par):
        '''
        Function where only sidebands are considered.

        Parameters
        ----------

        - x: function variable
        - par: function parameters
        '''
        if self.removePeak and TMath.Abs(x[0] - self.peakMass) < self.peakDelta:
            TF1.RejectPoint()
            return 0
        if self.removeSecPeak and TMath.Abs(x[0] - self.secPeakMass) < self.secPeakDelta:
            TF1.RejectPoint()
            return 0

        return self._ExpIntegralNorm(x, par)

    def GetBkgSideBandsFunc(self, integral):
        numPar = 2
        funcBkgSB = TF1('bkgSBfunc', self._SideBandsFunc, self.minMass, self.maxMass, numPar)
        funcBkgSB.SetParNames('BkgInt', 'Slope')
        funcBkgSB.SetParameters(integral, -2.)
        funcBkgSB.SetLineColor(kBlue+2)
        return funcBkgSB

    def GetBkgFullRangeFunc(self, integral):
        numPar = 2
        funcBkg = TF1('bkgFunc', self._ExpIntegralNorm, self.minMass, self.maxMass, numPar)
        funcBkg.SetParNames('BkgInt', 'Slope')
        funcBkg.SetParameters(integral, -2.)
        funcBkg.SetLineColor(kGreen+2)
        return funcBkg
