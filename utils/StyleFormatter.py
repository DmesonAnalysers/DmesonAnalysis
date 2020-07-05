'''
Script with helper methods for style settings
'''

from ROOT import gStyle, gROOT, TGaxis # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kWhite, kGray, kRed, kBlue, kGreen # pylint: disable=import-error,no-name-in-module
from ROOT import kMagenta, kAzure, kCyan, kOrange, kYellow, kSpring, kTeal, kViolet, kPink # pylint: disable=import-error,no-name-in-module
from ROOT import kFullCircle, kFullSquare, kFullDiamond, kFullCross, kFullTriangleUp, kFullTriangleDown # pylint: disable=import-error,no-name-in-module
from ROOT import kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kOpenTriangleUp, kOpenTriangleDown # pylint: disable=import-error,no-name-in-module

# pylint: disable=too-many-branches, too-many-statements
def SetGlobalStyle(**kwargs):
    '''
    Method to set global style.

    Parameters
    ----------

    - padrightmargin (float), default = 0.035
    - padleftmargin (float), default = 0.12
    - padtopmargin (float), default = 0.035
    - padbottommargin (float), default = 0.12

    - titlesize (float), default = 0.050
    - titlesizex (float), default = 0.050
    - titlesizey (float), default = 0.050
    - titlesizez (float), default = 0.050

    - labelsize (float), default = 0.045
    - labelsizex (float), default = 0.045
    - labelsizey (float), default = 0.045
    - labelsizez (float), default = 0.045

    - titleoffset (float), default = 1.2
    - titleoffsetx (float), default = 1.2
    - titleoffsey (float), default = 1.2
    - titleoffsetz (float), default = 1.2

    - opttitle (int), default = 0
    - optstat (int), default = 0

    - padtickx (int), default = 1
    - padticky (int), default = 1

    - maxdigits (int), default no max value

    - palette (int), default kBird
    '''

    # pad margins
    if 'padrightmargin' in kwargs:
        gStyle.SetPadRightMargin(kwargs['padrightmargin'])
    else:
        gStyle.SetPadRightMargin(0.035)

    if 'padleftmargin' in kwargs:
        gStyle.SetPadLeftMargin(kwargs['padleftmargin'])
    else:
        gStyle.SetPadLeftMargin(0.12)

    if 'padtopmargin' in kwargs:
        gStyle.SetPadTopMargin(kwargs['padtopmargin'])
    else:
        gStyle.SetPadTopMargin(0.035)

    if 'padbottommargin' in kwargs:
        gStyle.SetPadBottomMargin(kwargs['padbottommargin'])
    else:
        gStyle.SetPadBottomMargin(0.1)

    # title sizes
    if 'titlesize' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesize'], 'xyz')
    else:
        gStyle.SetTitleSize(0.050, 'xyz')

    if 'titlesizex' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesizex'], 'x')
    if 'titlesizey' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesizex'], 'y')
    if 'titlesizez' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesizex'], 'z')

    # label sizes
    if 'labelsize' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsize'], 'xyz')
    else:
        gStyle.SetLabelSize(0.045, 'xyz')

    if 'labelsizex' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsizex'], 'x')
    if 'labelsizey' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsizex'], 'y')
    if 'labelsizez' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsizex'], 'z')

    # title offsets
    if 'titleoffset' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffset'], 'xyz')
    else:
        gStyle.SetTitleOffset(1.2, 'xyz')

    if 'titleoffsetx' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffsetx'], 'x')
    if 'titleoffsety' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffsety'], 'y')
    if 'titleoffsetz' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffsetz'], 'z')

    # other options
    if 'opttitle' in kwargs:
        gStyle.SetOptTitle(kwargs['opttitle'])
    else:
        gStyle.SetOptTitle(0)

    if 'optstat' in kwargs:
        gStyle.SetOptStat(kwargs['optstat'])
    else:
        gStyle.SetOptStat(0)

    if 'padtickx' in kwargs:
        gStyle.SetPadTickX(kwargs['padtickx'])
    else:
        gStyle.SetPadTickX(1)

    if 'padticky' in kwargs:
        gStyle.SetPadTickY(kwargs['padticky'])
    else:
        gStyle.SetPadTickY(1)

    gStyle.SetLegendBorderSize(0)

    if 'maxdigits' in kwargs:
        TGaxis.SetMaxDigits(kwargs['maxdigits'])

    if 'palette' in kwargs:
        gStyle.SetPalette(kwargs['palette'])

    gROOT.ForceStyle()


def SetObjectStyle(obj, **kwargs):
    '''
    Method to set root object style.

    Parameters
    ----------

    - obj: object to set style

    - linecolor (int) default 1 (black)
    - linealpha (float) default 1
    - linewitdh (int) default 2
    - linestyle (int) default 1

    - markercolor (int) default 1 (black)
    - markeralpha (float) default 1
    - markerstyle (int) default 20 (full circle)
    - markersize (int) default 20 (full circle)

    - fillcolor (int) default no filling
    - fillalpha (float) default 1
    - fillstyle (int) default 0 (no style)

    - color (int) sets same color for line, marker and fill
    - alpha (float) sets same alpha for line, marker and fill
    '''

    # alpha parameters
    lalpha = kwargs.get('linealpha', 1)
    malpha = kwargs.get('markeralpha', 1)
    falpha = kwargs.get('fillalpha', 1)
    if 'alpha' in kwargs:
        lalpha = kwargs['alpha']
        malpha = kwargs['alpha']
        falpha = kwargs['alpha']
    if 'linealpha' in kwargs:
        lalpha = kwargs['linealpha']
    if 'markeralpha' in kwargs:
        malpha = kwargs['markeralpha']
    if 'fillalpha' in kwargs:
        falpha = kwargs['fillalpha']

    # line styles
    if 'linecolor' in kwargs:
        if lalpha < 1:
            obj.SetLineColorAlpha(kwargs['linecolor'], lalpha)
        else:
            obj.SetLineColor(kwargs['linecolor'])
    else:
        if lalpha < 1:
            obj.SetLineColorAlpha(1, lalpha)
        else:
            obj.SetLineColor(1)

    if 'linewidth' in kwargs:
        obj.SetLineWidth(kwargs['linewidth'])
    else:
        obj.SetLineWidth(2)

    if 'linestyle' in kwargs:
        obj.SetLineStyle(kwargs['linestyle'])
    else:
        obj.SetLineStyle(1)

    # marker styles
    if 'markercolor' in kwargs:
        if malpha < 1:
            obj.SetMarkerColorAlpha(kwargs['markercolor'], malpha)
        else:
            obj.SetMarkerColor(kwargs['markercolor'])
    else:
        if malpha < 1:
            obj.SetMarkerColorAlpha(1, malpha)
        else:
            obj.SetMarkerColor(1)

    if 'markersize' in kwargs:
        obj.SetMarkerSize(kwargs['markersize'])
    else:
        obj.SetMarkerSize(1)

    if 'markerstyle' in kwargs:
        obj.SetMarkerStyle(kwargs['markerstyle'])
    else:
        obj.SetMarkerStyle(20)

    # fill styles
    if 'fillcolor' in kwargs:
        if falpha < 1:
            obj.SetFillColorAlpha(kwargs['fillcolor'], falpha)
        else:
            obj.SetFillColor(kwargs['fillcolor'])

    if 'fillstyle' in kwargs:
        obj.SetFillStyle(kwargs['fillstyle'])

    #global color
    if 'color' in kwargs:
        if lalpha < 1:
            obj.SetLineColorAlpha(kwargs['color'], lalpha)
        else:
            obj.SetLineColor(kwargs['color'])
        if malpha < 1:
            obj.SetMarkerColorAlpha(kwargs['color'], malpha)
        else:
            obj.SetMarkerColor(kwargs['color'])
        if falpha < 1:
            obj.SetFillColorAlpha(kwargs['color'], falpha)
        else:
            obj.SetFillColor(kwargs['color'])


def DivideCanvas(canv, nPads):
    '''
    Method to divide ROOT canvases

    Parameters
    ----------

    - canv: TCanvas to be divided
    - nPads: number (int) of pads in which divide the canvas

    '''
    if nPads < 2:
        canv.cd()
    elif nPads in [2, 3]:
        canv.Divide(int(nPads), 1)
    elif nPads in [4, 6, 8]:
        canv.Divide(int(nPads/2), 2)
    elif nPads in [5, 7]:
        canv.Divide(int((nPads+1)/2), 2)
    elif nPads in [12, 15]:
        canv.Divide(int(nPads/3), 3)
    elif nPads in [10, 11]:
        canv.Divide(4, 3)
    elif nPads in [13, 14]:
        canv.Divide(5, 3)
    elif 15 < nPads <= 20 and nPads % 4 == 0:
        canv.Divide(int(nPads/4), 4)
    elif 15 < nPads <= 20 and nPads % 4 != 0:
        canv.Divide(5, 4)
    elif nPads == 21:
        canv.Divide(7, 3)
    elif 21 < nPads <= 25:
        canv.Divide(5, 5)
    elif nPads > 25 and nPads % 2 == 0:
        canv.Divide(int(nPads/2), 2)
    else:
        canv.Divide(int((nPads+1)/2), 2)


def GetROOTColor(color='kBlack'):
    '''
    Method to retrieve a ROOT color

    Parameters
    ----------

    - color: color according to ROOT TColor convention

    Returns
    ----------

    - ROOT color corresponding to input color

    '''
    cMapROOT = {'kBlack': kBlack, 'kWhite': kWhite, 'kGrey': kGray,
                'kRed': kRed, 'kBlue': kBlue, 'kGreen': kGreen,
                'kTeal': kTeal, 'kAzure': kAzure, 'kCyan': kCyan,
                'kOrange': kOrange, 'kYellow': kYellow, 'kSpring': kSpring,
                'kMagenta': kMagenta, 'kViolet': kViolet, 'kPink': kPink}

    ROOTcolor = None
    for colorKey in cMapROOT:
        if colorKey in color:
            ROOTcolor = cMapROOT.get(colorKey)
            break
    if ROOTcolor:
        for shade in range(0, 11):
            if f' + {shade}' in color or f'+{shade}' in color:
                ROOTcolor += shade
                break
            elif f' - {shade}' in color or f'-{shade}' in color:
                ROOTcolor -= shade
                break

    return ROOTcolor


def GetROOTMarker(marker='kFullCircle'):
    '''
    Method to retrieve the ROOT marker map

    Parameters
    ----------

    - color: color according to ROOT TColor convention

    Returns
    ----------

    - ROOT color corresponding to input color

    '''
    mMapROOT = {'kFullCircle': kFullCircle, 'kFullSquare': kFullSquare, 'kFullDiamond': kFullDiamond,
                'kFullCross': kFullCross, 'kFullTriangleUp': kFullTriangleUp, 'kFullTriangleDown': kFullTriangleDown,
                'kOpenCircle': kOpenCircle, 'kOpenSquare': kOpenSquare, 'kOpenDiamond': kOpenDiamond,
                'kOpenCross': kOpenCross, 'kOpenTriangleUp': kOpenTriangleUp, 'kOpenTriangleDown': kOpenTriangleDown}

    if marker in mMapROOT:
        ROOTmarker = mMapROOT.get(marker)
    else:
        ROOTmarker = None

    return ROOTmarker
