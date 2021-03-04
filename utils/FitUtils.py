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
class BkgFitFuncCreator:
    '''
    Class to handle custom background functions as done by AliHFInvMassFitter. Mainly designed
    to provide functions for sidebands fitting

    Parameters
    -------------------------------------------------
    - funcName: function to use. Currently implemented: 'expo', 'pol0', 'pol1', 'pol2', 'pol3'
    - minMass:  lower extreme of fitting interval
    - maxMass:  higher extreme of fitting interval
    - numSigmaSideBands: number of widths excluded around the peak
    - peakMass: peak mass (if not defined the signal region will not be excluded from the function)
    - peakSigma: peak width
    - secPeakMass: second peak mass (if not defined the second-peak region will not be excluded from the function)
    - secPeakSigma: second peak width
    '''
    __implFunc = {'expo': '_ExpoIntegralNorm',
                  'pol0': '_Pol0IntegralNorm',
                  'pol1': '_Pol1IntegralNorm',
                  'pol2': '_Pol2IntegralNorm',
                  'pol3': '_Pol3IntegralNorm'
                  }

    __numPar = {'expo': 2,
                'pol0': 1,
                'pol1': 2,
                'pol2': 3,
                'pol3': 4
                }

    def __init__(self, funcName, minMass, maxMass, numSigmaSideBands=0., peakMass=0.,
                 peakSigma=0., secPeakMass=0., secPeakSigma=0.):
        if funcName not in self.__implFunc:
            raise ValueError(f'Function \'{funcName}\' not implemented')
        self.funcName = funcName
        self.minMass = minMass
        self.maxMass = maxMass
        self.peakMass = peakMass
        self.peakDelta = peakSigma * numSigmaSideBands
        self.secPeakMass = secPeakMass
        self.secPeakDelta = secPeakSigma * numSigmaSideBands
        self.funcSBCallable = None
        self.funcFullCallable = None

        self.removePeak = False
        self.removeSecPeak = False
        if self.peakMass > 0. and self.peakDelta > 0.:
            self.removePeak = True
        if self.secPeakMass > 0. and self.secPeakDelta > 0.:
            self.removeSecPeak = True

    def _ExpoIntegralNorm(self, x, par):
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

    def _Pol0IntegralNorm(self, x, par): # pylint: disable=unused-argument
        '''
        Constant Function normalized to its integral.

        Parameters
        ----------
        - x: function variable
        - par: function parameters
            par[0]: normalisation (integral of background)
        '''
        return par[0] / (self.maxMass - self.minMass)

    def _Pol1IntegralNorm(self, x, par):
        '''
        Linear function normalized to its integral.
        See AliHFInvMassFitter::FitFunction4Bkg for more information.

        Parameters
        ----------
        - x: function variable
        - par: function parameters
            par[0]: normalisation (integral of background)
            par[1]: angular coefficient
        '''
        return par[0] / (self.maxMass - self.minMass) + par[1] * (x[0] - 0.5 * (self.maxMass + self.minMass))

    def _Pol2IntegralNorm(self, x, par):
        '''
        Second order polinomial function normalized to its integral.
        See AliHFInvMassFitter::FitFunction4Bkg for more information.

        Parameters
        ----------
        - x: function variable
        - par: function parameters
            par[0]: normalisation (integral of background)
            par[1]: a
            par[2]: b
        '''
        firstTerm = par[0] / (self.maxMass - self.minMass)
        secondTerm = par[1] * (x[0] - 0.5 * (self.maxMass + self.minMass))
        thirdTerm = par[2] * (x[0]**2 - 1 / 3. * (self.maxMass**3 - self.minMass**3) / (self.maxMass - self.minMass))
        return firstTerm + secondTerm + thirdTerm

    def _Pol3IntegralNorm(self, x, par):
        '''
        Third order polinomial function normalized to its integral.
        See AliHFInvMassFitter::FitFunction4Bkg for more information.

        Parameters
        ----------
        - x: function variable
        - par: function parameters
            par[0]: normalisation (integral of background)
            par[1]: a
            par[2]: b
            par[3]: c
        '''
        firstTerm = par[0] / (self.maxMass - self.minMass)
        secondTerm = par[1] * (x[0] - 0.5 * (self.maxMass + self.minMass))
        thirdTerm = par[2] * (x[0]**2 - 1 / 3. * (self.maxMass**3 - self.minMass**3) / (self.maxMass - self.minMass))
        fourthTerm = par[3] * (x[0]**3 - 1 / 4. * (self.maxMass**4 - self.minMass**4) / (self.maxMass - self.minMass))
        return firstTerm + secondTerm + thirdTerm + fourthTerm

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

        return getattr(self, self.__implFunc[self.funcName])(x, par)

    def GetSideBandsFunc(self, integral):
        '''
        Return the ROOT.TF1 function defined on the sidebands

        Parameters
        --------------------------------------
        integral: integral of the histogram to fit, obtained with TH1.Integral('width')

        Returns
        ---------------------------------------
        funcBkgSB: ROOT.TF1
            Background function
        '''
        self.funcSBCallable = self._SideBandsFunc # trick to keep away the garbage collector
        funcBkgSB = TF1('bkgSBfunc', self.funcSBCallable, self.minMass, self.maxMass, self.__numPar[self.funcName])
        funcBkgSB.SetParName(0, 'BkgInt')
        funcBkgSB.SetParameter(0, integral)
        for iPar in range(1, self.__numPar[self.funcName]):
            funcBkgSB.SetParameter(iPar, 1.)
        funcBkgSB.SetLineColor(kBlue+2)
        return funcBkgSB

    def GetFullRangeFunc(self, func):
        '''
        Return the ROOT.TF1 function defined on the full range

        Parameters
        --------------------------------------
        func: function from GetSideBandsFunc() after the histogram fit

        Returns
        ---------------------------------------
        funcBkg: ROOT.TF1
            Background function
        '''
        self.funcFullCallable = getattr(self, self.__implFunc[self.funcName]) # trick to keep away the garbage collector
        funcBkg = TF1('bkgFunc', self.funcFullCallable, self.minMass, self.maxMass, self.__numPar[self.funcName])
        funcBkg.SetParName(0, 'BkgInt')
        for iPar in range(0, self.__numPar[self.funcName]):
            funcBkg.SetParameter(iPar, func.GetParameter(iPar))
        funcBkg.SetLineColor(kGreen+2)
        return funcBkg
