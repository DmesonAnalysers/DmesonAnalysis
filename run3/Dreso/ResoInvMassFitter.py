import ROOT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mplhep
import sys
import argparse
import yaml
import os
from pathlib import Path


def fit_invariant_mass_with_background(histogram, x_min, x_max):
    """
    Fit a Gaussian model with a polynomial background to an invariant mass histogram using ROOFit.
    
    Parameters:
    - histogram: ROOT.TH1, histogram to be fitted.
    - x_min: float, minimum x-value (mass) to consider in the fit.
    - x_max: float, maximum x-value (mass) to consider in the fit.
    
    Returns:
    - workspace: ROOT.RooWorkspace, workspace containing the model and fit result.
    - mass_frame: ROOT.RooPlot, frame containing the fit and histogram for plotting.
    """
    # Create RooRealVar for the mass
    mass = ROOT.RooRealVar("mass", "Invariant Mass", x_min, x_max)

    # Convert TH1 histogram to a ROOT.RooDataHist
    data = ROOT.RooDataHist("data", "Dataset with invariant mass", ROOT.RooArgList(mass), histogram)

    # Define a Gaussian model for the signal
    mean = ROOT.RooRealVar("mean", "Mean of Gaussian", 2.53, 2.52, 2.56) #Ds1
    # mean = ROOT.RooRealVar("mean", "Mean of Gaussian", 0.524, 0.51, 0.54) #Ds1_delta
    # mean = ROOT.RooRealVar("mean", "Mean of Gaussian", 0.145, 0.142, 0.148) #Dstar
    sigma = ROOT.RooRealVar("sigma", "Width of Gaussian", 0.01, 0., 0.20)
    width = ROOT.RooRealVar("width", "Width of BW", 0.00046) #Ds1
    # width = ROOT.RooRealVar("width", "Width of BW", 0.0000415) #Dstar

    voigt = ROOT.RooVoigtian("voigt", "Voigtian PDF", mass, mean, width, sigma)

    #Define Generig PDF for Background
    mthr = ROOT.RooRealVar("mthr", "Threshold mass", 2.5076) #Ds1_
    # mthr = ROOT.RooRealVar("mthr", "Threshold mass", 0.498) #Ds1_delta
    # mthr = ROOT.RooRealVar("mthr", "Threshold mass", 0.1385) # Dstar
    l = ROOT.RooRealVar("l", "expo coefficient", 5, 0, 1000)
    fbkg = ROOT.RooGenericPdf("fbkg", "resonance bkg PDF", "sqrt(abs(mass-mthr))*exp(-l*(mass-mthr))", ROOT.RooArgSet(mass, mthr, l))

    # Combine signal and background into a composite model
    signal_fraction = ROOT.RooRealVar("signal_fraction", "Fraction of Signal", 0.01, 0, 0.5)
    Nsig = ROOT.RooRealVar("Nsig","Number of signal events",0.01*histogram.Integral(), 0 ,0.3*histogram.Integral())
    Nbkg = ROOT.RooRealVar("Nbkg","Number of backgrouund events",0.9*histogram.Integral(), 0 ,histogram.Integral())
    total_pdf = ROOT.RooAddPdf("total_pdf", "Total PDF with Nsig as parameter", ROOT.RooArgList(voigt, fbkg), ROOT.RooArgList(Nsig, Nbkg))

    # Fit the composite model to the data
    fit_result = total_pdf.fitTo(data, ROOT.RooFit.Save())

    # Create a frame to draw the fit result and data
    mass_frame = mass.frame()
    data.plotOn(mass_frame)
    total_pdf.plotOn(mass_frame)
    total_pdf.plotOn(mass_frame, ROOT.RooFit.Components("fbkg"), ROOT.RooFit.LineStyle(ROOT.kDashed))

    # Save to workspace
    workspace = ROOT.RooWorkspace(f"w__{histogram.GetName()}", f"workspace_{histogram.GetName()}")
    getattr(workspace, 'import')(total_pdf)
    getattr(workspace, 'import')(data)
    getattr(workspace, 'import')(fit_result, "fitResults")  # Give a name to your fit results for later access

    # Draw the frame
    canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    mass_frame.Draw()
    canvas.SaveAs("uglyplot2.pdf")
    return workspace, mass_frame

def roofit_plot_with_matplotlib(axis, histogram, workspace, x_min, x_max, fit_min, fit_max, plot):
    """
    Plot a TH1 and ROOFit fitted PDFs using matplotlib .
    
    Parameters:
    - axis: plt.axis, axis where to draw the plot.
    - histogram: ROOT.TH1, Invariant mass histogram
    - x_min: float, minimum of mass plot range.
    - x_max: float, maximum of mass plot range.
    - fit_min: float, minimum of mass fit range (for normalization porpuses).
    - fit_max: float, maximum of mass fit range (for normalization porpuses).
    
    Returns:
    - workspace: ROOT.RooWorkspace, workspace containing the model and fit result.
    - mass_frame: ROOT.RooPlot, frame containing the fit and histogram for plotting.
    """

    # Access mass variable and data
    mass = workspace.var("mass")
    data = workspace.data("data")
    # frac = workspace.var("signal_fraction")
    nsig = workspace.var("Nsig")
    nbkg = workspace.var("Nbkg")
    fit_result = workspace.obj("fitResults")
   
    # Create the underlying histogram
    Nmin = histogram.FindBin(x_min)
    Nmax = histogram.FindBin(x_max)
    N = Nmax - Nmin
    N_fitmin = histogram.FindBin(fit_min)
    N_fitmax = histogram.FindBin(fit_max)
    N_fit = N_fitmax - N_fitmin
    N_plot = 250
    bin_center = np.zeros(N)
    bin_content = np.zeros(N)
    bin_error = np.zeros(N)
    # Accessing each data point
    for i in range(N):
        bin_center[i] = histogram.GetBinCenter(Nmin + i)
        bin_content[i] = histogram.GetBinContent(Nmin + i)
        bin_error[i] = histogram.GetBinError(Nmin + i)

    totPDF = workspace.pdf("total_pdf")
    signalPDF = workspace.pdf("voigt")
    bkgPDF = workspace.pdf("fbkg")
    y_totPDF = np.zeros(N_fit)
    y_signalPDF = np.zeros(N_fit)
    y_bkgPDF = np.zeros(N_fit)
    x_fit = np.linspace(fit_min, fit_max, N_fit)
    x_plot = np.linspace(fit_min, fit_max, N_plot)
    y_tot_plot = np.zeros(N_plot)
    y_signal_plot = np.zeros(N_plot)
    y_bkg_plot = np.zeros(N_plot)
    for i, x in enumerate(x_fit):
        mass.setVal(x)
        y_totPDF[i] = totPDF.getVal(ROOT.RooArgSet(mass)) 
        y_signalPDF[i] = signalPDF.getVal(ROOT.RooArgSet(mass)) 
        y_bkgPDF[i] = bkgPDF.getVal(ROOT.RooArgSet(mass)) 
    for i, x in enumerate(x_plot):
        mass.setVal(x)
        y_tot_plot[i] = totPDF.getVal(ROOT.RooArgSet(mass)) 
        y_signal_plot[i] = signalPDF.getVal(ROOT.RooArgSet(mass)) 
        y_bkg_plot[i] = bkgPDF.getVal(ROOT.RooArgSet(mass)) 

    # Retrieving useful information for legend 
    chi2 = totPDF.createChi2(data)
    mean = workspace.var("mean")
    sigma = workspace.var("sigma")
    width = workspace.var("width")
    
    # Compute S, B and Significance
    hwhm = 0.5346 * width.getVal() + (0.2166 * width.getVal()**2 + 2*np.log(2)*sigma.getVal()**2)**0.5 #from JQSRT 17, P233
    mass.setRange("intRange", mean.getVal() - 3*hwhm , mean.getVal() + 3*hwhm)
    bkg_3HWHM = bkgPDF.createIntegral(ROOT.RooArgSet(mass), ROOT.RooFit.NormSet(ROOT.RooArgSet(mass)), ROOT.RooFit.Range("intRange"))
    bkg = nbkg.getVal()*bkg_3HWHM.getVal()
    bkg_err = nbkg.getError()*bkg_3HWHM.getVal()
    sig = nsig.getVal()
    sig_err = nsig.getError()
    sig_plus_bkg = nsig.getVal()+nbkg.getVal()*bkg_3HWHM.getVal() 
    if sig <= 0 or bkg <= 0:
        return [0, 0, 0, 0]
    if chi2.getVal() < 2* (N_fit-fit_result.floatParsFinal().getSize()):
        significance = nsig.getVal()/np.sqrt(sig_plus_bkg)
        significance_err = significance*np.sqrt((sig_err**2 + bkg_err**2) / (4. * sig_plus_bkg**2) + (bkg/sig_plus_bkg) * sig_err**2 / sig**2)
    else:
        significance = 0
        significance_err = 0
        sig = 0
        sig_err = 0
    if plot:
        #Create legend text
        textstr1 = '\n'.join(( 
            '\u03BC = %.1f \u00B1 %.1f MeV/c$^2$' % (mean.getVal([0])*1000, mean.getError()*1000),
            '\u03C3 =  %.1f \u00B1 %.1f MeV/c$^2$' % (sigma.getVal([0])*1000, sigma.getError()*1000),
            '\u0393 =  %.1f \u00B1 %.1f MeV/c$^2$' % (2*width.getVal([0])*1000, width.getError()*1000),
            'HWHM = %.1f MeV/c$^2$'  % (hwhm * 1000),
            'S = %d \u00B1 %d' % (int(nsig.getVal()), int(nsig.getError())),
            'B(3 HWHM) = %d \u00B1 %d' % (int(nbkg.getVal()*bkg_3HWHM.getVal()), int((nbkg.getError()*bkg_3HWHM.getVal()))),
            'B = %d \u00B1 %d' % (int(nbkg.getVal()), int(nbkg.getError())),
                ))
        textstr2 = '\n'.join(( 
            '',
            '\u03C7$^2$/ndf = %.1f/%d' % (chi2.getVal(), N_fit-fit_result.floatParsFinal().getSize()),
            'Significance(3 HWHM) = %.1f \u00B1 %.1f' % (significance, significance_err),
                ))
        axis.text(0.015, 0.95, textstr2, transform=axis.transAxes, fontsize=15, verticalalignment='top')

        # Plotting data and fit PDFs
        axis.errorbar(bin_center, bin_content, yerr=bin_error, fmt='o', label='Data', markersize=5, color='black')
        
        norm = histogram.Integral(N_fitmin,N_fitmax)/np.sum(y_totPDF)
        f = nsig.getVal()/((nbkg.getVal() + nsig.getVal()))  
        
        # if pdg_reso == 10433:
        #     name_reso = "Ds1plus"
        #     label_mass = r"M($\mathrm{D}^{*+}\mathrm{K_S^0}$)$"
        #     pdg_v0 = 310
        #     pdg_d = 413
        #     extra_info_loc = ["lower right", "upper center"]
        # elif pdg_reso == 435:
        #     name_reso = "Ds2starplus"
        #     label_mass = r"M($\mathrm{D}^+\mathrm{K_S^0}$)$"
        #     pdg_v0 = 310
        #     pdg_d = 411
        #     extra_info_loc = ["upper left", "lower right"]
        # label_mass = r"M($\mathrm{D}^+\mathrm{K_S^0}$)$ GeV/{c^2}"

        axis.plot(x_plot, (1-f)*norm*y_bkg_plot, label='Bkg', linewidth=2, linestyle='--', color='red')
        axis.plot(x_plot, f*norm*y_signal_plot, label='Sig', linewidth=1, color='blue')
        axis.plot(x_plot, norm*y_tot_plot, label='Sig + Bkg', linewidth=2, color='blue')
        axis.fill_between(x_plot, f*norm*y_signal_plot, color='skyblue', alpha=0.4)
        axis.legend(fontsize=15)

        # Plot Options
        axis.set_title(f"{histogram.GetName()}", fontsize=18)
        # axis.set_xlabel(rf"{label_mass} (GeV/$c^2$)", fontsize=18)
        axis.set_xlabel(f"Invariant mass (GeV/c²)", fontsize=18)
        axis.set_ylabel(f'Counts/{((bin_center[2]-bin_center[1])*1000):.1f} MeV/c²', fontsize=18)
        axis.set_ylim(0., max(bin_content) * 1.5)
        props = dict(boxstyle='square', facecolor='white', alpha=0.8, edgecolor='black', zorder=10)
        axis.text(0.75, 0.4, textstr1, transform=axis.transAxes, fontsize=15, verticalalignment='top', color='black', bbox=props, backgroundcolor='white')

    results = [int(sig), int(sig_err), significance , significance_err ]
    return results
  
# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
    parser.add_argument('--output_directory','-d',default='.',
                    help='Directory for output plots. If ".", save in the same directory of the input file.')
    args = parser.parse_args()
    
    ROOT.gROOT.SetBatch(True)  # Run in batch mode; don't display canvases

    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    
    # # Configure matplotlib to use LaTeX
    # matplotlib.rcParams['text.usetex'] = True
    # matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    
    infilename = inputCfg['inputfile']
    outfilename = inputCfg['outputfile']
    dirname = inputCfg['dirname']
    histname = inputCfg['histname']

    input_file = ROOT.TFile(infilename, "read")
    folder = input_file.GetDirectory(dirname)
    th2_invMass = folder.Get(histname)

    outfile = ROOT.TFile(outfilename, "RECREATE")
    
    mplhep.style.use("ATLAS")
    rows = 2 
    cols = 2
    fig, ax = plt.subplots(nrows= rows, ncols= cols, figsize = (26,14))
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.3, wspace=0.3)

    for iPt, (ptMin, ptMax) in enumerate(zip(inputCfg['pt_bins']['min'], inputCfg['pt_bins']['max'])):

        proj = th2_invMass.ProjectionY(f"{th2_invMass.GetName()}_{inputCfg['pt_bins']['min'][iPt]}-{inputCfg['pt_bins']['max'][iPt]} GeV/c", th2_invMass.GetXaxis().FindBin(ptMin + 1e-6* ptMin), th2_invMass.GetXaxis().FindBin(ptMax - 1e-6* ptMax))
        outfile.cd()
        proj.Write()
        proj.Rebin(inputCfg['rebin'][iPt])
        workspace, mass_frame = fit_invariant_mass_with_background(proj, inputCfg["fit_range"]["min"], inputCfg["fit_range"]["max"])
        irow = int(iPt / cols)
        icol = iPt % cols
        results = roofit_plot_with_matplotlib(ax[irow, icol], proj, workspace, inputCfg["mass_range"]["min"], inputCfg["mass_range"]["max"], inputCfg["fit_range"]["min"], inputCfg["fit_range"]["max"], plot = True)
    outfile.Close()
    fig.savefig(f"{args.output_directory}/{inputCfg['plotname']}.png")

