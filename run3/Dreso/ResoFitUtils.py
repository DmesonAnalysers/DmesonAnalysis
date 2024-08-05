'''
utility for fitting with Roofit D mesons particle spectra
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ROOT
from particle import Particle

def fit_invariant_mass_with_background(histogram, x_min, x_max, pdg_id, delta_m):
    """
    Fit invariant mass histogram using ROOFit.
    
    Parameters:
    - histogram: ROOT.TH1, histogram to be fitted.
    - x_min: float, minimum x-value (mass) to consider in the fit.
    - x_max: float, maximum x-value (mass) to consider in the fit.
    
    Returns:
    - workspace: ROOT.RooWorkspace, workspace containing the model and fit result.
    - mass_frame: ROOT.RooPlot, frame containing the fit and histogram for plotting.
    """
    supported_ids = [10433, 435, 411] #[Ds1, Ds2*, D+]

    if pdg_id not in supported_ids:
        raise KeyError(f"pdg_id {pdg_id} not supported!")
    
    mass = ROOT.RooRealVar("mass", "Invariant Mass", x_min, x_max)
    data = ROOT.RooDataHist("data", "Dataset with invariant mass", ROOT.RooArgList(mass), histogram)


    if pdg_id == 10433:
        # voigtian signal and exp-pow bkg
        m_r = Particle.from_pdgid(pdg_id).mass*1e-3
        w_r = Particle.from_pdgid(pdg_id).width/2*1e-3
        m_d = Particle.from_pdgid(413).mass*1e-3
        m_v0 = Particle.from_pdgid(310).mass*1e-3
        if delta_m:
            mean = ROOT.RooRealVar("mean", "Mean of Gaussian", m_r - m_d, m_r - m_d - 0.02, m_r - m_d + 0.02)
            mthr = ROOT.RooRealVar("mthr", "Threshold mass", m_v0)
        else:
            mean = ROOT.RooRealVar("mean", "Mean of Gaussian", m_r, m_r - 0.02, m_r + 0.02)
            mthr = ROOT.RooRealVar("mthr", "Threshold mass", m_d + m_v0)
        sigma = ROOT.RooRealVar("sigma", "Width of Gaussian", 0.01, 0., 0.02)
        width = ROOT.RooRealVar("width", "Width of BW", w_r)
        signal_pdf = ROOT.RooVoigtian("signal_pdf", "Voigtian PDF", mass, mean, width, sigma)
        l = ROOT.RooRealVar("l", "expo coefficient", 5, 0, 1000)
        bkg_pdf = ROOT.RooGenericPdf("bkg_pdf", "resonance bkg PDF", "sqrt(abs(mass-mthr))*exp(-l*(mass-mthr))", ROOT.RooArgSet(mass, mthr, l))
    
    elif pdg_id == 435:
        # voigtian signal and pol1 bkg
        m_r = Particle.from_pdgid(pdg_id).mass*1e-3
        w_r = Particle.from_pdgid(pdg_id).width/2*1e-3
        m_d = Particle.from_pdgid(411).mass*1e-3
        m_v0 = Particle.from_pdgid(310).mass*1e-3
        if delta_m:
            mean = ROOT.RooRealVar("mean", "Mean of Gaussian", m_r - m_d)
        else:
            mean = ROOT.RooRealVar("mean", "Mean of Gaussian", m_r, m_r - 0.01, m_r + 0.01)
        sigma = ROOT.RooRealVar("sigma", "Width of Gaussian", 0.01, 0., 0.02)
        width = ROOT.RooRealVar("width", "Width of BW", w_r)
        signal_pdf = ROOT.RooVoigtian("signal_pdf", "Voigtian PDF", mass, mean, width, sigma)
        a0 = ROOT.RooRealVar("a0", "Constant Coefficient", 100, 0, 1e7)
        a1 = ROOT.RooRealVar("a1", "Linear Coefficient", 0, -1e3, 1e3)
        bkg_pdf = ROOT.RooPolynomial("bkg_pdf", "Polynomial Background", mass, ROOT.RooArgList(a0, a1))

    elif pdg_id == 411:
        # Gaussian signal and pol3 bkg
        m = Particle.from_pdgid(411).mass*1e-3
        mean = ROOT.RooRealVar("mean", "Mean of Gaussian", m, m - 0.02, m + 0.02)
        sigma = ROOT.RooRealVar("sigma", "Width of Gaussian", 0.01, 0., 0.02)
        signal_pdf = ROOT.RooGaussian("signal_pdf", "Gaussian Model", mass, mean, sigma)
        a0 = ROOT.RooRealVar("a0", "Constant Coefficient", 100, 0, 1e7)
        a1 = ROOT.RooRealVar("a1", "Linear Coefficient", 0, -1e3, 1e3)
        a2 = ROOT.RooRealVar("a2", "Quadratic Coefficient", 0, -1e3, 1e3)
        a3 = ROOT.RooRealVar("a3", "Cubic Coefficient", 0, -1e3, 1e3)
        bkg_pdf = ROOT.RooPolynomial("bkg_pdf", "Polynomial Background", mass, ROOT.RooArgList(a0, a1, a2, a3))
        width = ROOT.RooRealVar("width", "Width of BW", 0)

    # Combine signal and background into a composite model
    Nsig = ROOT.RooRealVar("Nsig","Number of signal events",0.01*histogram.Integral(), 0 ,0.3*histogram.Integral())
    Nbkg = ROOT.RooRealVar("Nbkg","Number of backgrouund events",0.9*histogram.Integral(), 0 ,histogram.Integral())
    total_pdf = ROOT.RooAddPdf("total_pdf", "Total PDF with Nsig as parameter", ROOT.RooArgList(signal_pdf, bkg_pdf), ROOT.RooArgList(Nsig, Nbkg))

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
    getattr(workspace, 'import')(width)
    getattr(workspace, 'import')(fit_result, "fitResults")  # Give a name to your fit results for later access

    # Draw the frame
    canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    mass_frame.Draw()
    canvas.SaveAs("uglyplot.pdf")
    return workspace, mass_frame

def roofit_plot_with_matplotlib(axis, histogram, workspace, x_min, x_max, fit_min, fit_max, plot):
    """
    Plot a TH1 and ROOFit fitted PDFs using matplotlib .
    
    Parameters:
    - axis: plt.axis, axis where to draw the plot.
    - histogram: ROOT.TH1, Invariant mass histogramselectedMasses = invMassReso[(ptReso>ptMin) & (ptReso<ptMax)]
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
    signalPDF = workspace.pdf("signal_pdf")
    bkgPDF = workspace.pdf("bkg_pdf")
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
        # sig = 0
        # sig_err = 0
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
        axis.text(0.75, 0.3, textstr1, transform=axis.transAxes, fontsize=15, verticalalignment='top', color='black', bbox=props, backgroundcolor='white')

    results = [int(sig), int(sig_err), significance , significance_err ]
    return results